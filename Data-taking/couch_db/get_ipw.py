# need to source the tellie environment first

from orca_side import tellie_database

# set up access to database
database = tellie_database.TellieDatabase('http://couch.snopl.us', 'telliedb', 'tellieadmin', 'IsThisPINReal42')

# load all docs for tuning view
channelStart = 1
channelEnd = 2
startKey = [channelStart, {}]
endKey = [channelEnd, {}]
rows = database.db.view('_design/tuning/_view/tuning', startkey=startKey, endkey=endKey, include_docs=True)

for row in rows:
	row_data = row.doc
	print row

