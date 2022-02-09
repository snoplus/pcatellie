# need to source the tellie environment first

from orca_side import tellie_database
import json

# set up access to database
database = tellie_database.TellieDatabase('http://couch.snopl.us', 'telliedb', 'tellieadmin', 'IsThisPINReal42')

#loop through files and upload

for i in range (1,96):

	# load file
	file_name = "all_docs/"+str(i)+".json"
	with open(file_name) as json_file:
	    data = json.load(json_file)

	# upload file		
	database.save(data)
	
	print "Done: ", i

print "DONE!"

