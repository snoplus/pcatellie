function(doc){
    if(doc.type=="TUNING"){
        emit([doc.channel, doc.ipw, doc.run_range, doc.version, doc.pass], [doc.timestamp]);
    }
}
