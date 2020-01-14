
// Raw data conversion for the BDAQ53 producer, 2018
// jorge.duarte.campderros@cern.ch

#define DEBUG_1 0

#include "eudaq/RD53ADecoder.hh"

#include "eudaq/DataConverterPlugin.hh"
#include "eudaq/StandardEvent.hh"
#include "eudaq/Utils.hh"
#include "eudaq/RawDataEvent.hh"
#include "eudaq/Timer.hh"
#include "eudaq/Logger.hh"

#if USE_LCIO
#  include "IMPL/LCEventImpl.h"
#  include "IMPL/TrackerRawDataImpl.h"
#  include "IMPL/LCCollectionVec.h"
#  include "lcio.h"
#  include "IMPL/TrackerDataImpl.h"
#  include "IMPL/LCGenericObjectImpl.h"
#  include "UTIL/CellIDEncoder.h"
#endif

#if USE_EUTELESCOPE
// Eudaq
#  include "EUTelRD53ADetector.h"
#  include "EUTelRunHeaderImpl.h"
// Eutelescope
#  include "EUTELESCOPE.h"
#  include "EUTelSetupDescription.h"
#  include "EUTelEventImpl.h"
#  include "EUTelTrackerDataInterfacerImpl.h"
#  include "EUTelGenericSparsePixel.h"
using eutelescope::EUTELESCOPE;
#endif

#include <stdlib.h>
#include <cstdint>
#include <map>
#include <vector>
#include <array>
#include <memory>

#ifdef DEBUG_1
#include <bitset>
#include <iomanip>
#endif 

namespace eudaq 
{
  // The event type for which this converter plugin will be registered
  // Modify this to match your actual event type (from the Producer)

  //static const char* EVENT_TYPE = "BDAQ53"; // Jun 2018
  static const char* EVENT_TYPE   = "bdaq53a"; // Oct 2018

#if USE_LCIO && USE_EUTELESCOPE
  static const int chip_id_offset = 30;
#endif

  // Declare a new class that inherits from DataConverterPlugin
  class BDAQ53ConverterPlugin : public DataConverterPlugin
  {
  public:

    bool ldb = 0;

    // This is called once at the beginning of each run.
    // You may extract information from the BORE and/or configuration
    // and store it in member variables to use during the decoding later.
    virtual void Initialize(const Event & bore, const Configuration & cnf) 
    {
#ifndef WIN32
      (void)cnf; // just to suppress a warning about unused parameter cnf
#endif
      // XXX: Is there any parameters? 
      _get_bore_parameters(bore);
      //bore.Print(std::cout);
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // This should return the trigger ID (as provided by the TLU)
    // if it was read out, otherwise it can either return (unsigned)-1,
    // or be left undefined as there is already a default version.
    virtual unsigned GetTriggerID(const Event & ev) const 
    {
      // Make sure the event is of class RawDataEvent
      if(const RawDataEvent * rev = dynamic_cast<const RawDataEvent *> (&ev)) {
	if( rev->IsBORE() || rev->IsEORE() ) {
	  return 0;
	}
	else if( rev->NumBlocks() < 1 ) {
	  return (unsigned)-1;
	}
	else if(rev->GetBlock(0).size() == 0) {
	  return (unsigned)-1;
	}

	if( ldb )
	  std::cout << "BDAQ53ConverterPlugin GetStandardSubEvent rev->NumBlocks() " << rev->NumBlocks()
		    << ", rev->GetBlock(0).size " << rev->GetBlock(0).size()
		    << std::endl << std::flush;

	unsigned int trigger_number = RD53ADecoder::get_trigger_number( rev->GetBlock(0) );

	if( trigger_number == static_cast<unsigned>(-1) ) {
	  EUDAQ_ERROR("No trigger word found in the first 32-block.");
	  return (unsigned)-1;
	}
	else {
	  if( ldb )
	    std::cout << "trigger_number " << trigger_number  << std::endl << std::flush;
	  return trigger_number;
	}
      }
      // Not proper event or not enable to extract Trigger ID
      return (unsigned)-1;
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Here, the data from the RawDataEvent is extracted into a StandardEvent.
    // The return value indicates whether the conversion was successful.

    virtual bool GetStandardSubEvent( StandardEvent & sev, const Event & ev ) const 
    {
      if( ev.IsBORE() || ev.IsEORE() ) {
	// nothing to do
	return true;
      }

      // Clear memory
      _data_map->clear();
                
      // Number of event headers: must be, at maximum, the number of
      // subtriggers of the Trigger Table (TT)
      //unsigned int n_event_headers = 0;
                
      // If we get here it must be a data event
      const RawDataEvent & ev_raw = dynamic_cast<const RawDataEvent &>(ev);

      // Trigger number is the same for all blocks (sensors)
      uint32_t const trigger_number = GetTriggerID(ev);

      if( static_cast<uint32_t>(trigger_number) == -1 ) {

	// Dummy plane

	StandardPlane plane( 0, EVENT_TYPE, "RD53A" );

	plane.SetSizeZS( RD53A_NCOLS, RD53A_NROWS, 0, RD53A_MAX_TRG_ID,
			StandardPlane::FLAG_DIFFCOORDS | StandardPlane::FLAG_ACCUMULATE );

	sev.AddPlane( plane );

	return false;
      }

      //else if(trigger_number != 0)
      //{
      //    // Got a trigger word
      //    event_status |= E_EXT_TRG;
      //}

      // [JDC] each block could corresponds to one sensor attach to 
      //       the card. So far, we are using only 1-sensor, but it 
      //       could be extended at some point

      if( ldb )
	std::cout << "BDAQ53ConverterPlugin::GetStandardSubEvent: ev_raw.NumBlocks() = " << ev_raw.NumBlocks()
		  << std::endl << std::flush;

      for( size_t i = 0; i < ev_raw.NumBlocks(); ++i ) {

	// Create a standard plane representing one sensor plane
	StandardPlane plane( i, EVENT_TYPE, "RD53A" );

	if( ldb )
	  std::cout << "  got plane " << i << std::endl <<std::flush;

	// ZS = zero suppressed data
	// Each frame can have different coordinates
	// accumulate multiple frames

	plane.SetSizeZS( RD53A_NCOLS, RD53A_NROWS, 0, RD53A_MAX_TRG_ID,
			StandardPlane::FLAG_DIFFCOORDS | StandardPlane::FLAG_ACCUMULATE );

	const auto & decoded_data = get_decoded_data( ev_raw.GetBlock(i), i );

	if( ldb )
	  std::cout << "  decoded_data.n_event_headers() " << decoded_data.n_event_headers()
		    << std::endl << std::flush;

	plane.SetTLUEvent( trigger_number );

	// fill data-header (be used?)

	if( decoded_data.n_event_headers() < 4294967295 ) { // 2^32-1

	  for( unsigned int dh = 0; dh < decoded_data.n_event_headers(); ++dh ) {

	    uint32_t trig_id = decoded_data.trig_id()[dh];
	    //decoded_data->header_number;
#if DEBUG_1
	    uint32_t bcid = decoded_data.bcid()[dh];
	    uint32_t trig_tag = decoded_data.trig_tag()[dh];
	    std::cout << "EH  " << std::setw(9) << bcid
		      << std::setw(9) << trig_id << std::setw(9) 
		      << trig_tag << " " << std::bitset<32>(data_word) << std::endl;
#endif
	    if( ldb )
	      std::cout << "  frame " << dh
			<< ", trig " << trig_id
			<< ", BC " << decoded_data.bcid()[dh]
		;

	    for( const auto & hits: decoded_data.hits(dh) ) {
	      // pixel = {col, row, ToT, pivot, frame }
	      plane.PushPixel( hits[0], hits[1], hits[2], false, trig_id );
	      if( ldb )
		std::cout << " : " << hits[0] << "." << hits[1] << "." << hits[2];
	    }

	    if( ldb )
	      std::cout << std::endl <<std::flush;

	  } // dh

	}

	sev.AddPlane( plane );

      } // blocks (planes)

      return true;

    }

#if USE_LCIO && USE_EUTELESCOPE
    // This is where the conversion to LCIO is done
    virtual lcio::LCEvent * GetLCIOEvent(const Event * /*ev*/) const 
    {
      return 0;
    }
            
    virtual bool GetLCIOSubEvent(lcio::LCEvent & lcioEvent, const Event & eudaqEvent) const
    {
      //std::cout << "getlciosubevent (I4) event " << eudaqEvent.GetEventNumber() << " | " << GetTriggerID(eudaqEvent) << std::endl;
      if(eudaqEvent.IsBORE() ) {
	// shouldn't happen
	return true;
      } 
      else if( eudaqEvent.IsEORE() ) {
	// nothing to do
	return true;
      }

      // Clear the memory
      _data_map->clear();
                
      // set type of the resulting lcio event
      lcioEvent.parameters().setValue(eutelescope::EUTELESCOPE::EVENTTYPE, eutelescope::kDE);
      // pointer to collection which will store data
      LCCollectionVec * zsDataCollection = nullptr;
      // it can be already in event or has to be created
      bool zsDataCollectionExists = false;
      try 
	{
	  zsDataCollection = static_cast<LCCollectionVec*>(lcioEvent.getCollection("zsdata_rd53a"));
	  zsDataCollectionExists = true;
	} 
      catch(lcio::DataNotAvailableException& e) 
	{
	  zsDataCollection = new LCCollectionVec(lcio::LCIO::TRACKERDATA);
	}
      //	create cell encoders to set sensorID and pixel type
      CellIDEncoder<TrackerDataImpl> zsDataEncoder(eutelescope::EUTELESCOPE::ZSDATADEFAULTENCODING, zsDataCollection);
                
      // this is an event as we sent from Producer
      // needs to be converted to concrete type RawDataEvent
      const RawDataEvent & ev_raw = dynamic_cast <const RawDataEvent &>(eudaqEvent);
      std::vector<eutelescope::EUTelSetupDescription*> setupDescription;
                
      // [JDC - XXX: Just in case we're able to readout more chips
      for(size_t chip = 0; chip < ev_raw.NumBlocks(); ++chip ) {
	//const std::vector <unsigned char> & buffer=dynamic_cast<const std::vector<unsigned char> &>(ev_raw.GetBlock(chip));
	const RawDataEvent::data_t & raw_data = ev_raw.GetBlock(chip);
	if (lcioEvent.getEventNumber() == 0) {
	  // XXX - Hardcoded pitch (50x50 um)
	  eutelescope::EUTelPixelDetector * currentDetector = new eutelescope::EUTelRD53ADetector(0.05,0.05);
	  currentDetector->setMode("ZS");
	  setupDescription.push_back( new eutelescope::EUTelSetupDescription(currentDetector)) ;
	}
	zsDataEncoder["sensorID"] = ev_raw.GetID(chip)+chip_id_offset;  
	zsDataEncoder["sparsePixelType"] = eutelescope::kEUTelGenericSparsePixel;

	// prepare a new TrackerData object for the ZS data
	// it contains all the hits for a particular sensor in one event
	std::unique_ptr<lcio::TrackerDataImpl > zsFrame(new lcio::TrackerDataImpl);
	// set some values of "Cells" for this object
	zsDataEncoder.setCellID(zsFrame.get());

	// this is the structure that will host the sparse pixel
	// it helps to decode (and later to decode) parameters of all hits (x, y, charge, ...) to
	// a single TrackerData object (zsFrame) that will correspond to a single sensor in one event
	std::unique_ptr< eutelescope::EUTelTrackerDataInterfacerImpl<eutelescope::EUTelGenericSparsePixel> >
	  sparseFrame(new eutelescope::EUTelTrackerDataInterfacerImpl<eutelescope::EUTelGenericSparsePixel>(zsFrame.get()));

	// Reassemble full 32-bit FE data word: 
	// The word is constructed using two 32b+32b words (High-Low) 
	// In order to build the FE high and low words, you need
	// 4 blocks of uint8 (unsigned char) on little endian for the FE high word
	// and 4 blocks of uint8 more (on little endian) for the FE low word
	// (see _reassemble_word)
	const auto & decoded_data = get_decoded_data(ev_raw.GetBlock(chip),chip);
	for(unsigned int dh=0; dh < decoded_data.n_event_headers(); ++dh ) {
	  for(const auto & hits: decoded_data.hits(dh) ) {
	    sparseFrame->emplace_back(hits[0],hits[0],hits[2],decoded_data.trig_id()[dh]);
	  }
	}
	// write TrackerData object that contains info from one sensor to LCIO collection
	zsDataCollection->push_back(zsFrame.release());
      }
      // add this collection to lcio event
      if( (!zsDataCollectionExists )  && (zsDataCollection->size() != 0) ) {
	lcioEvent.addCollection( zsDataCollection, "zsdata_rd53a" );
      }
      if (lcioEvent.getEventNumber() == 0 ) {
	// do this only in the first event
	LCCollectionVec * apixSetupCollection = nullptr;
	bool apixSetupExists = false;
	try 
	  {
	    apixSetupCollection=static_cast<LCCollectionVec*>(lcioEvent.getCollection("rd53a_setup"));
	    apixSetupExists = true;
	  }
	catch(...) 
	  {
	    apixSetupCollection = new LCCollectionVec(lcio::LCIO::LCGENERICOBJECT);
	  }

	for(auto const & plane: setupDescription ) {
	  apixSetupCollection->push_back(plane);
	}
	if(!apixSetupExists) {
	  lcioEvent.addCollection(apixSetupCollection, "rd53a_setup");
	}
      }
      return true;
    }
#endif

  private:

    void _get_bore_parameters(const Event & )
    {
      // XXX Just if need to initialize paremeters from the BORE
      //std::cout << bore << std::endl;
    }

    //--------------------------------------------------------------------------
    const RD53ADecoder & get_decoded_data( const RawDataEvent::data_t & raw_data,const int & chip) const
    {
      // memorize if not there
      if(_data_map->find(chip) == _data_map->end() ) {
	_data_map->emplace(chip,RD53ADecoder(raw_data));
      }
      return _data_map->at(chip); //_data_map->find(chip)->second;
    }

    //--------------------------------------------------------------------------
    // The decoded data memory (pointer to avoid constantness of the Get** functions)
    std::unique_ptr<std::map<int,RD53ADecoder> >  _data_map;

    //--------------------------------------------------------------------------
    // The constructor can be private, only one static instance is created
    // The DataConverterPlugin constructor must be passed the event type
    // in order to register this converter for the corresponding conversions
    // Member variables should also be initialized to default values here.
    BDAQ53ConverterPlugin() : 
      DataConverterPlugin(EVENT_TYPE),
      _data_map(std::unique_ptr<std::map<int,RD53ADecoder>>(new std::map<int,RD53ADecoder>)) 
    {
    }

    // The single instance of this converter plugin
    static BDAQ53ConverterPlugin m_instance;

  };

  // Instantiate the converter plugin instance
  BDAQ53ConverterPlugin BDAQ53ConverterPlugin::m_instance;

} // namespace eudaq
