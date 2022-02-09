import os
import sys

def load_env():
    ratdb_user = os.getenv("RATDB_USER")
    ratdb_pw = os.getenv("RATDB_PW")
    ratdb_add = os.getenv("RATDB_ADD")
    tables_loc = os.getenv("TABLES_LOC")
    return ratdb_user, ratdb_pw, ratdb_add, tables_loc

if __name__=="__main__":

    if len(sys.argv) != 2:
        print "Please provide the radtb file name..."
        exit()
    else:
        ratdb_file_name = sys.argv[1]
        if str(ratdb_file_name[-6:]) != ".ratdb":
            print "Please provide a valid ratdb file..."
            exit()
    print ratdb_file_name

    # load env
    ratdb_user, ratdb_pw, ratdb_add, tables_loc = load_env()

    # prepare command
    com = "ratdb -s postgres://" + str(ratdb_user) + ":" + str(ratdb_pw) + "@" + str(ratdb_add) + " upload " + str(tables_loc + ratdb_file_name)
    #os.system(com)
