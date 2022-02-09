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
        
    def read_pin_sequence(self):
        print "Call to read_pin_sequence(), returning random mean and rms values"
        if (np.random.rand(1)[0] < 0.5):
            return int(np.random.rand(1)[0]*65355), float(np.random.rand(1)[0]*250), self._chan
        return
                   
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
