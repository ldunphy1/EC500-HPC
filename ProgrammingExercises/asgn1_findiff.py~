import math
def func(x):
	return math.sin(x)
def derivative(x):
	return math.cos(x)
def forward_diff(x, h):
	return (func(x+h)-func(x))/h
def backward_diff(x, h):
	return (func(x)-func(x-h))/h
def central_diff(x, h):
	return (func(x+h/2)-func(x-h/2))/h
h = 1.0
x = 1.0
#Add: set the output precision to 15 decimal places (the default is 6)!
#Loop over 17 values of h: 1 to 10^(-16).
for i in range(0, 17):
	h *= 0.1
	#Print h, the forward, backward, and central difference of sin(x) at 1,
	#as well as the exact derivative cos(x) at 1.
	#Should you use spaces or tabs as delimiters? Does it matter?
	print '%.15e %.15f %.15f %.15f %.15f\n' % (h, forward_diff(x,h), backward_diff(x,h), central_diff(x,h), math.cos(1.0))
