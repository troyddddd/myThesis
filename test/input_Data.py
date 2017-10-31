


with open ("test_data.txt",'r') as f:
	for line in f:
		temp1,temp2,temp3 = (ele.strip() for ele in line.split(" ")[:3])
	print temp1,temp2,temp3


