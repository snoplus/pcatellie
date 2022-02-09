import ROOT
import json
from rat import ratdb

def load_fibre_data(fibres):

    Us = []
    Vs = []
    Ws = []
    fib_names = []

    ROOT.RAT.DB.Get().LoadDefaults()

    for fibre in fibres:
        entry = ROOT.RAT.DB.Get().GetLink("FIBRE", fibre)
        vector = ROOT.TVector3(entry.GetD("u"), entry.GetD("v"), entry.GetD("w"))
        #print "Direction is ( %.10f | %.10f | %.10f ) mm" % (vector.X(), vector.Y(), vector.Z())
        fib_names.append( fibre )
        Us.append( vector.X() )
        Vs.append( vector.Y() )
        Ws.append( vector.Z() )
    return fib_names, Us, Vs, Ws

if __name__=="__main__":

    fibres = ["FT001A", "FT002A", "FT003A", "FT004A", "FT005A", "FT006A", "FT007A", "FT008A", "FT009A", "FT010A", "FT011A", "FT012A", "FT013A", "FT014A", "FT015A", "FT016A", "FT017A", "FT018A", "FT019A", "FT020A", "FT021A", "FT022A", "FT023A", "FT024A", "FT025A", "FT026A", "FT027A", "FT028A", "FT029A", "FT030A", "FT031A", "FT032A", "FT033A", "FT034A", "FT035A", "FT036A", "FT037A", "FT038A", "FT039A", "FT040A", "FT041A", "FT042A", "FT043A", "FT044A", "FT045A", "FT046A", "FT047B", "FT048A", "FT049A", "FT050A", "FT051A", "FT052A", "FT053A", "FT054A", "FT055A", "FT056A", "FT057A", "FT058A", "FT059A", "FT060A", "FT061A", "FT062A", "FT063A", "FT064A", "FT065A", "FT066A", "FT067A", "FT068A", "FT069A", "FT070A", "FT071A", "FT072A", "FT073A", "FT074A", "FT075A", "FT076A", "FT077A", "FT078A", "FT079A", "FT080A", "FT081A", "FT082A", "FT083A", "FT084A", "FT085A", "FT086A", "FT087A", "FT088A", "FT089A", "FT090A", "FT091A", "FT092A", "FT093A", "FT094A", "FT101A"]

    fib_names, U, V, W = load_fibre_data(fibres)

    print fib_names
    print U
    print V
    print W

    json_data = {
    "fibre": fib_names,
    "u": U,
    "v": V,
    "w": W
    }

    # Save to ratdb file
    with open("default_directions.ratdb", "w") as f:
        json.dump(json_data, f)
