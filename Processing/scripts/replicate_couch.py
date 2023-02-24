import couchdb
import os

def load_env():
    c_add = os.getenv("COUCH_ADD")
    c_db = os.getenv("COUCH_DB")
    c_user = os.getenv("COUCHDB_USER")
    c_pw = os.getenv("COUCHDB_PW")
    return c_add, c_db, c_user, c_pw

def get_all_docs():
	couchserver = couchdb.Server("http://"+c_user+":"+c_pw+"@"+c_add)
	db = couchserver[c_db]

	rows = db.view('_all_docs', include_docs=True)
	data = [row['doc'] for row in rows]

	return data

def upload_to_chdb(data):
    couchserver_target = couchdb.Server("http://snoplusdb:*****@couch.snopl.us") ### put pw here
    db_target = couchserver_target[c_db]

    counter = 0
    for source_doc in data:
    	try:
        	# save the new doc
        	del source_doc['_id']
        	del source_doc['_rev']
        	#print source_doc
        	doc = db_target.save(source_doc)

        	counter += 1
        	if counter%10 == 0: 
        		print "Done " + str(counter) + " docs!"

    	except:
        	print "Couldn't create the doc :/"
        	exit()

if __name__=="__main__":
	c_add, c_db, c_user, c_pw = load_env()
	print c_add, c_db, c_user, c_pw 
	data = get_all_docs()

	print "Total docs: ", len(data)

	upload_to_chdb(data)