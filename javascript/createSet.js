  // use an Object to create a unique array;
  function createSet(ai) {		
        var n = {}, bi=[];		
        for(var i=0; i< ai.length; i++) {		
            if (!n[ai[i]]) {		
                n[ai[i]] = true; 		
                bi.push(ai[i]); 		
            }		
        }		
        return bi;		
    }
