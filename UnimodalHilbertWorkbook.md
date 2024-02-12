# The Weak Lefschetz property and unimodality of Hilbert functions of random monomial algebras 
This page contains Macaulay2 code for generating random monomial algebras and the study of unimodality of Hilbert functions and whether or not the Weak Lefschetz property (WLP) holds. 


## Prerequisite 

We rely on the Macaulay2 package for generating Erdos-Renyi-type random monomial ideals. 
```
loadPackage "RandomMonomialIdeals"
```

The following small function that checks how many times a sign change occurred in a sequence of numbers. 
```
signChanges = (L) ->(
    c := 0;
    lastsign := 1;
    s := 1;
    for j from 1 to #L-1 do(
    	if L#j-L#(j-1) == 0 then s = s else s = (L#j-L#(j-1))/abs(L#j-L#(j-1));
	if lastsign == -s then(
	    c = c + 1; 
	    lastsign = s);
       	);
    return c;
    )
```

The next function checks if the WLP holds: 
```
wlp = (I) ->( 
    J = I + ideal sum flatten entries vars ring I;
    B = (ring J)^1/J;
    A = (ring I)^1/I;

    hB:={};
    i:=0;
    while hilbertFunction(i,B)=!=0 do(
	hB = append(hB,hilbertFunction(i,B));
	i=i+1;
	);
    hB = append(hB,0);
     
    hA:={};
    scan(#hB,i -> hA = append(hA,hilbertFunction(i,A)));

    -- check if hB_j == hA(j)-hA(j-1), the WLP condition: 
    --     apply(1..#hB-1, j->  hB_j == hA_j-hA_(j-1)) 
    firstIndexOfFailureWLP = scan(1..#hB-2, j-> if  hB_j =!= hA_j-hA_(j-1) then break j);  -- saving it in case we prefer to return this index. 
    -- if WLP holds SO FAR until the last index, then firstIndexOfFailureWLP will be null; therefore:
    -- now check at j= #hB-1, which is the first value where h_B(j)=0: 
    if firstIndexOfFailureWLP===null then firstIndexOfFailureWLP = hA_(#hB-1) <= hA_(#hB-2) else firstIndexOfFailureWLP = false;

    -- if WLP holds then firstIndexOfFailureWLP will be not null: 
    return firstIndexOfFailureWLP
    )
``` 


## Generating random monomial algebras  and testing unimodality of Hilbert function and the WLP 

Here is some Macaulay2 code we used to generate a few examples.
```
n=5;
R = QQ[x_1..x_n]
D=50;
apply(1..D,numer -> (
  p = toRR numer/D;
  print "p="; 
  print p; 

  N=100; 
  myData =  randomMonomialIdeals(n,D,p,N); 
  -- -*for verification*- dimStats(myData, ShowTally=>true);

  zeroDimSample = {};
  apply(myData, i->
    if (dim i)==0 then zeroDimSample = append(zeroDimSample,sub(i,R))
    );
  -- -*for verification*- dimStats zeroDimSample; 

  --here are all the h-vectors: "hvalues" is a HashTable, keys are ideals, values are h;
  --    and "multiplicities" is another one, same keys, values are number of times we saw that ideal in the simulation: 
  multiplicities = new MutableHashTable; 

  hvalues = new MutableHashTable; 
    for I in zeroDimSample do(
	  -- check if multiplicities#I exists, if it does do not recompute h!, but increase multiplicty. 
	  if multiplicities #? I  then ( 
	    --record every time we see this same ideal I: 
	    multiplicities#I = multiplicities#I+1; 
	    ) else ( 
	    --print "didn't find I in keys"; 
	    h := {};
	    for i from 0 to D-1 do (	    
		h = append(h,hilbertFunction(i,I))
		);	
	    hvalues# I=h;
	    multiplicities#I=1; -- count this at least once - create the key
	    );

	);
  hvalues = new HashTable from hvalues;
  multiplicities = new HashTable from multiplicities;
  -* 
    -- And here is the non-hash-table version of the same computation:
    hvalues = {};
      for I in zeroDimSample do(
  	h := {};
  	for i from 0 to D-1 do (
  	    h = append(h,hilbertFunction(i,I)));
  	hvalues = append(hvalues,h);
  	);
  *- 

  --what proportion are unimodal?  answer: fractionUnimodal = countNonUnimodal/#zeroDimSample
  countNonUnimodal = 0;
  countNonUnimodalUnique = 0; 
  countWLP = 0; 
  countWLPunique = 0; 
  scanKeys(hvalues,I-> if signChanges hvalues#I > 1 then (
	  countNonUnimodal = countNonUnimodal+1*multiplicities#I;
  	countNonUnimodalUnique = countNonUnimodalUnique+1; 
	  ) else (  -- for the UNIMODAL ones, check if the ideal I has WLP! 
	  if wlp(I) then (
	    countWLP=countWLP+1*multiplicities#I;
	    countWLPunique = countWLPunique+1; 
	    ) 
	  )
    );
  -* 
    -- the non-hash-table version:
    apply(hvalues,h-> if signChanges h > 1 then (
  	  countNonUnimodal=countNonUnimodal+1;
    	)
    );
  *- 
```

Now that we have computed some samples, it would be good to write the data to a file: 
```
f = openOutAppend "unimodalStatistics.txt";
f << endl; 
f << countNonUnimodal << "|" << #zeroDimSample <<  "|" << toRR countNonUnimodal/#zeroDimSample;  
f << " -- for D" << D << " n" << n << " p" << p;
if countNonUnimodal > 0 then f << "; and of these, "<< countNonUnimodalUnique << " unique ideals with non-unimodal Hilbert functions.";
f << endl; 
close f;

f = openOutAppend "WLPStatistics.txt";
f << endl; 
f << countWLP << "|" << (#zeroDimSample-countNonUnimodal) <<  "|" << toRR countWLP/(#zeroDimSample-countNonUnimodal);  
f << " -- for D" << D << " n" << n << " p" << p; 
f << "; and of these, "<< countWLPunique << " unique ideals with WLP." << endl; 
close f;
```

### Expected Hilbert computations

We wish to compute $E[h_d]$ for $d=0,\dots,D+1$ for the zero-dimensional ideals we consider.
```
h = new MutableHashTable; 
apply(0..D-1, d-> 
    (
    h#d={};
    scanValues(hvalues, hi-> 
	h#d = append(h#d, hi_d) -- "hi" is one of the LISTS of h-vectors stored in the hashtable hvalues; hi_d is then the d-th entry, thus the value of h_d for some ideal in our sample.
    )
)
);
-* -- the non-hash-table version:
apply(0..D-1, d-> 
    (
    h#d={};
    apply(0..#hvalues-1, i-> 
	h#d = append(h#d, hvalues_i_d) -- hvalues_i is the ith h-vector LIST. then _d is the d-th entry of that. 
    )
)
);
*- 
```
The result  of the "apply/scanValues" loop above is that  $h$ is a hash table such that h#d is the list of all values of h_d we see in the data.

Now, $E[h_d]$ is simply the mean of this list:
```
expectedHilbert = new MutableHashTable; 
apply(0..D-1, d-> expectedHilbert #d = toRR(sum h#d / #h#d)); -- mean := (sum fData)/s.SampleSize;
-- at this point, "values expectedHilbert" is a list of numbers, the d-th entry is the expected value of h_d. 
-- WARNING: THE CODE ASSUMES THAT THE KEYS OF THE HASHTABLE expectedHilbert ARE SORTED IN INCREASING ORDER!!!! 
```

Let us check if the expected Hilbert function is unimodal:
```
signChanges values expectedHilbert 
```
and if we wish to write this to a file as well: 
```
f = openOutAppend "unimodalStatistics.txt";
f << "E[h] unimodal: " << if(signChanges values expectedHilbert)<=1 then "yes" else "no" << endl;
close f;
```

## Random Artinian algebras 

### Model 1: $I_{ER} + (x_1^D, ..., x_n^D)$

``` 
n=4
D = 50; 
N=50;  -- sample size 

-*   -- just going 0.1,0.2,...0.9:
apply(9,numer-> (
 p = toRR 0.1+0.1*numer; 
*- 
-- when ER has dim 1 .. n (not just zero): 
apply(1..n,numer -> (	
 p = toRR 1/(D^numer);
 
print "p="; 
print p; 

mySample = randomArtinAlgebras((n,D,p,N)); -- this generates hashtable with Hvectors and corresponding Ideals as keys! 

f = openOutAppend "unimodalStatistics.txt";
f << endl; 
f << "Number of non-unimodal h-vectors in the sample: " << (modeStats mySample#Hvectors)#Nonunimodal ;
f << " -- for D" << D << " n" << n << " p" << p << endl;  
f << "E[h] unimodal: " << if(signChanges (modeStats mySample#Hvectors)#Means )<=1 then "yes" else "no" << endl;
close f;


f = openOutAppend "WLPStatistics.txt";
f << endl; 
f << "rate of sample with wlp = " << (wlpStats mySample#Ideals)#Rate << " -- for D" << D << " n" << n << " p" << p; 
close f;

)
)
```


If you run the code sufficiently many times, you may notice a pattern. For example, when D= 100 (not when D=50), around p=1/D^2, it seems that there is a low probability of WLP holding. Let us test this on a sample: 

```
n=4
D=50
apply(1..3,numer -> (	
 D = numer * 50; 
 p = toRR 1/(D^2);
print "D=";
print D;  
print "p="; 
print p; 

N=50;
mySample = randomArtinAlgebras((n,D,p,N)); -- this generates hashtable with Hvectors and corresponding Ideals as keys! 

print "Number of non-unimodal h-vectors in the sample: "; 
print (modeStats mySample#Hvectors)#Nonunimodal ;
print "E[h] unimodal: ";
if(signChanges (modeStats mySample#Hvectors)#Means )<=1 then print "yes" else  print "no"; 

print "rate of sample with wlp = ";
print (wlpStats mySample#Ideals)#Rate ;

)
)
```

### Model 2: $I_{ER} + (x_1, ..., x_n)^D$

This is the exact same code as above, except the `mySample = randomArtinAlgebras((n,D,p,N));` line should be replaced by:
```
mySample = randomArtinPlusAlgebras((n,D,p,N)); 
```

### Random level algebras 

This is computing through the socle. 
We show one example.


``` 
n=3 
D = 100; 
N=100; 
p=0.9
N=20

mySample = randomLevelAlgebrasNontrivial3((D,p,N)); -- this generates hashtable with Hvectors and corresponding Ideals as keys! 

-- check unimodality of E[h]: 
signChanges (modeStats mySample#Hvectors)#Means 
-- if(signChanges (modeStats mySample#Hvectors)#Means )<=1 then "yes" else "no"  

print("% of sample with wlp = ");
(wlpStats mySample#Ideals)#Rate

```

