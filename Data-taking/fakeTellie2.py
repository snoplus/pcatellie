import time
import numpy as np
from SimpleXMLRPCServer import SimpleXMLRPCServer 
    
class tellie(object):

    def init_channel(self, chan, no_pulses, pulse_sep, trigger_delay, pulse_width, pulse_height, fibre_delay):
        self._chan = chan
        settings = {"channel" : chan,
                    "no_pulses" : no_pulses,
                    "pulse_delay" : pulse_sep,
                    "trigger_delay" : trigger_delay,
                    "pulse_width" : pulse_width,
                    "pulse_height" : pulse_height,
                    "fibre_delay" : fibre_delay}
        print "Called init_channel with %s" % (settings)
        return settings

    def set_pulse_number(self, no_pulses):
        print "Called set_pulse_number with: %i" % (int(no_pulses))
        return 0
    
    def trigger_averaged(self):
        time.sleep(5)
        print "Call to trigger_averaged()"
    
    def read_pin_sequence(self, ipw):
        def polyn(x):
            p0 = 535288
            p1 = -518.175
            p2 = 0.236674
            p3 = -5.56311e-05
            p4 = 6.87594e-09
            p5 = -4.25697e-13
            p6 = 1.04202e-17
            y = p6*pow(x,6) + p5*pow(x,5) + p4*pow(x,4) + p3*pow(x,3) + p2*pow(x,2) + p1*pow(x,1) + p0
            return y
        print "Call to read_pin_sequence(), returning mean and rms values"
        if (ipw < 3500):
            return 65535, 5.0, self._chan
        elif (ipw > 8000):
                return 0, 5.0, self._chan
        else:
            return polyn(ipw), 5.0, self._chan

    def fire_sequence(self):
        print "Call to fire_sequence()"
        
    def stop(self):
        print "Call to stop()" 

    def is_connected(self):
        return "Tellie server is responding"

    def test(self):
        return "Tellie server is responding"
                
server = SimpleXMLRPCServer(("localhost", 5030), logRequests=True, allow_none=True)

root = tellie()
server.register_instance(root, allow_dotted_names=True)

try:
    print 'Serving... Control-c to exit'
    server.serve_forever()
except KeyboardInterrupt:
    print 'Killing...'
    raise
