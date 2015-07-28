#!/usr/bin/env python 

# vary selection
s = [ 0, 0.05, 0.25, 0.5 ]

type = [0,1,2]

# vary size
nmito = [ 10, 20, 30, 40 ]
nmito_error = [ 1, 5, 10 ]

# ascus intercept
intercept = [ 0.05, 0.1, 0.25, 0.5 ]

nreplicates = 200;

exe = "/home/uccoaku/yeast_heteroplasmy/src/xyeast_unequal_ascus"

count = 0

# ascus slope
for s_i in s:
    for type_i in type:
        for nmito_i in nmito:
            for nmito_error_i in nmito_error:
                if nmito_error == nmito_i:
                    continue

                for intercept_i in intercept:
                    slope = 1.0 - 2 * intercept_i

                    for rep_i in range(0,nreplicates):

                        count += 1
                        print("echo " + str(count))

                        print(exe + " " + str(s_i) 
                                + " " + str(type_i) 
                                + " " + str(nmito_i) 
                                + " " + str(nmito_error_i) 
                                + " " + str(intercept_i) 
                                + " " + str(slope))
