import xmlrpclib
#from fakeTellie.py import tellie

s = xmlrpclib.ServerProxy('http://localhost:5030')
print s.is_connected()
ipw = 8100
settings = s.init_channel(1, 200000, 1, 900, ipw, 8000, 10)
print s.read_pin_sequence(settings['pulse_width'])
