import csv

channel = []
ipw = []

with open('setpoints.csv') as csvfile:
    data = list(csv.reader(csvfile))

channels = len(data)
for i in range (1, channels):
	channel.append(int(data[i][0]))
	ipw.append(int(data[i][1]))


	doc_text = '{ \n'
	doc_text += '\t "type": "TUNING",\n'
	doc_text += '\t "index": "",\n'
	doc_text += '\t "timestamp": "2019-08-028T04:45:52.40Z",\n'
	doc_text += '\t "version": 0,\n'
	doc_text += '\t "run_range": [1, 300000],\n'
	doc_text += '\t "pass": 0,\n'
	doc_text += '\t "channel": ' + data[i][0] + ',\n'
	doc_text += '\t "ipw": ' + data[i][1] + '\n'
	doc_text += '}'

	file_name = data[i][0]+".json"
	file = open(file_name, 'w')
	file.write(doc_text)
	file.close()

