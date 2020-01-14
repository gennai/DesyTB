#!/usr/bin/env python

import serial
import sys
import hameg as h
import time

def main():
    vbb_dst=float(sys.argv[1]) # value to which the back-bias should be ramped
    sour=serial.Serial(h.dev, h.baud_rate, rtscts=h.rtscts)
    ret_value = 0 # 0 = success, 1 = error

    for ivbb in range (0, len(h.vbb_chan)):
        vbb_src=h.meas_voltage(sour, h.vbb_chan[ivbb])
        if abs((vbb_dst - vbb_src)/(vbb_dst + 1.e-3)) < 1.e-3:
            continue
        h.ramp_from_to(sour, h.vbb_chan[ivbb], (vbb_src, vbb_dst), 100)
        time.sleep(1)
        vbb_meas=h.meas_voltage(sour, h.vbb_chan[ivbb])
        if abs((vbb_dst - vbb_meas)/(vbb_dst + 1.e-3)) > 1.e-3:
            ret_value = 1

    return ret_value


## execute the main
if __name__ == "__main__":
    sys.exit(main())
