

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


def read_pin_sequence(ipw):
    if (ipw < 3500):
        return 65535
    elif (ipw > 8000):
        return 0
    else:
        return (polyn(ipw))

print read_pin_sequence(1000)
print read_pin_sequence(4000)
print read_pin_sequence(5000)
print read_pin_sequence(6000)
print read_pin_sequence(7000)
print read_pin_sequence(8000)
print read_pin_sequence(9000)
print read_pin_sequence(10000)
print read_pin_sequence(11000)
