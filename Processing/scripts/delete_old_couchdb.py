import couchdb
import os

c_add = os.getenv("COUCH_ADD")
c_db = os.getenv("COUCH_DB")
c_user = os.getenv("COUCHDB_USER")
c_pw = os.getenv("COUCHDB_PW")

couch = couchdb.Server("http://"+c_user+":"+c_pw+"@"+c_add)
db = couch[c_db]

#view_str = 'val1'
#view_str = 'fits'
view_str = 'val2'

all = db.view('_design/' + view_str + '/_view/'+ view_str , descending=False)

IDs = []
for row in all:
    print row.key

    try:
        print "Deleting", db[row.id]
        #db.delete(db[row.id])
        print ""
    except:
		print "!!! Exception deleting data doc !!!"

print "DONE ALL"
