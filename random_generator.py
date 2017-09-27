
def generate_sequence(start,end,length):
	from random import randint
	temp = []
	for i in range(length):
		temp.append(randint(start,end))
	return temp