#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#Parameters 
1 - name of cipher
2 - sbox or prob (sbox_2, sbox_4)
3- Solve (GUROBI/CPLEX)
4 - If no addition of inequalities than "-" otherwise group in which need to add inequalities "10" or "50", if this is "0" then all inequalities will added with all
5 - Whether to add 2 inequality or 3 or 4 inequalties (it should not be 1), if parameter 4 is - then there is no need of parameter 5

python Ineq_Reduction.py ESLICE sbox GUROBI -
python Ineq_Reduction.py ESLICE prob GUROBI -
python Ineq_Reduction.py ESLICE sbox GUROBI 0 2 
python Ineq_Reduction.py ESLICE sbox_2 GUROBI -
python Ineq_Reduction.py ESLICE prob CPLEX 10 2
"""


import sys
import time
import math
import numpy as np
from numpy import array, hstack, ones, vstack
import cdd  #pip install pycddlib

IN_COLAB = False
GUROBI_EXISTS = False
CPLEX_EXISTS = False

try:
  import google.colab
  IN_COLAB = True
except:
  IN_COLAB = False

if (sys.argv[3] == "GUROBI"):
    try: 
    	import gurobipy
    	from gurobipy import GRB
    	GUROBI_EXISTS = True
    except:
    	GUROBI_EXISTS = False

if (sys.argv[3] == "CPLEX"):
	try: 
		from docplex.mp.model import Model 
		CPLEX_EXISTS = True
	except:
		CPLEX_EXISTS = False

#############if manually want to choose the solve then set these parameters###########################
# GUROBI_EXISTS = False
# import gurobipy
# from gurobipy import GRB
# CPLEX_EXISTS = True		
# from docplex.mp.model import Model 
##############################################################

if (GUROBI_EXISTS == False) and (CPLEX_EXISTS == False):
	sys.exit()

if (IN_COLAB==True):
  # Setup the Gurobi environment with the WLS license
  e = gurobipy.Env(empty=True)
  e.start()

file_name = "_"+sys.argv[1]+"_"+sys.argv[2]+"_"+sys.argv[3]

def compute_polytope_halfspaces(vertices):
    V = vstack(vertices)
    t = ones((V.shape[0], 1))  # first column is 1 for vertices
    tV = hstack([t, V])
    mat = cdd.Matrix(tV, number_type='fraction')
    mat.rep_type = cdd.RepType.GENERATOR
    P = cdd.Polyhedron(mat)
    bA = array(P.get_inequalities())
    for i in range(0,len(bA)):
        mask = np.array([float(val).is_integer() for val in bA[i]])
        mul_factor = 2
        while(mask.sum()!=V.shape[1]+1):
            bA[i] = mul_factor*bA[i]
            mul_factor += 1
            mask = np.array([float(val).is_integer() for val in bA[i]])
    if bA.shape == (0,):  # bA == []
        return bA
    # the polyhedron is given by b + A x >= 0 where bA = [b|A]
    b, A = array(bA[:, 0]), array(bA[:, 1:])
    
    return (A, b)


def MILP_Solve(ineq_list,impossible_diff_arr):
        start_3 = time.process_time()
		
        if (GUROBI_EXISTS==True):
            if (IN_COLAB==True):
              print("environment set")
              m = gurobipy.Model(env=e)
            else:
              m = gurobipy.Model()
            for i in range(0,len(ineq_list)):
                m.addVar(vtype=GRB.BINARY,name="z%s" % str(i))
            m.update()
            variables = m.getVars()
            for i in range(0,len(impossible_diff_arr)):
                if(i%1000==0):
                    print(str(i)+" points covered")
                ineq_solve_count = (np.multiply(np.array(impossible_diff_arr[i]),np.array(ineq_list))).sum(1)
                less_than_zero = np.where(ineq_solve_count<0)[0]
                m.addConstr(sum([variables[x] for x in less_than_zero]) >=1)
            m.setObjective(sum(m.getVars()),GRB.MINIMIZE)
            m.update()
            end_3 = time.process_time()
            print(">>Time>> TTime Taken to write lp model: " + str(end_3 - start_3))
            m.write("problem"+ file_name +".lp")
            m.optimize()
            m.write("problem"+ file_name +".sol")
            m.setParam(GRB.Param.SolutionNumber, 1)
            sol_arr = np.array(m.X)
            final_sol = np.where(sol_arr == 1)[0].tolist()
        
        if (CPLEX_EXISTS==True):
            m = Model()
            m.binary_var
            for i in range(0,len(ineq_list)):
                m.binary_var(name="z%s" % str(i))
            variables = [i for i in m.iter_binary_vars()]
            for i in range(0,len(impossible_diff_arr)):
                ineq_solve_count = (np.multiply(np.array(impossible_diff_arr[i]),np.array(ineq_list))).sum(1)
                less_than_zero = np.where(ineq_solve_count<0)[0]
                m.add_constraint(sum([variables[x] for x in less_than_zero]) >=1)
            m.set_objective("min",sum(variables))
            end_3 = time.process_time()
            print(">>Time>> TTime Taken to write lp model: " + str(end_3 - start_3))
            m.export("CPLEX_Problem"+file_name +".lp")
            sol = m.solve(log_output=True)
            sol.export("CPLEX_Solution" +file_name+".sol")
            print(sol.get_objective_value())
            sol_arr = np.array(sol.get_values(variables)) 
            final_sol = np.where(sol_arr == 1)[0].tolist()
        
        
        f = open("ineq"+ file_name +".txt","w")
        ineq_list = ineq_list.astype(int)
        count = 1
        ineq_str = ""
        for ineqlitity in final_sol:
            inequality_rotated = ineq_list[ineqlitity].tolist()[1:] + [ineq_list[ineqlitity].tolist()[0]] 
            f.write(str(inequality_rotated).replace("[","").replace("]",",\n"))    
            ineq_str = ineq_str + str(inequality_rotated).replace("[","").replace("]",",\n")
            
            print_str = str(count)+"&$\\;\\;\\;"
            count += 1
            for i in range(0,len(inequality_rotated)-1):
                if (inequality_rotated[i]>=0):
                    print_str += "+"
                if( i <= 3):
                    print_str += str(inequality_rotated[i]) + "*x_" + str(3-i) + "  "
                if( i >= 4 and  i <= 7):
                    print_str += str(inequality_rotated[i]) + "*y_" + str(3-(i%4)) + "  "
                if( i >= 8):
                    print_str += str(inequality_rotated[i]) + "*p_" + str(1-(i%2)) + "  "
            print_str +=  " \\geq -" + str(str(inequality_rotated[-1])) + "$ \\\\"
            print(print_str)
        f.close()
        print("\n"+ineq_str)     
        ineq_list_mod = []
        print(final_sol)
        for x in final_sol:
            ineq_list_mod.append(ineq_list[x])
        return ineq_list
        
    
def print_DDT(table):
    for row in range(len(table)):
        for col in range(len(table[row])):
            print(table[row][col],end='');
            if col == len(table[row])-1:
                print("\n");

if (sys.argv[1] == "LICID"):
  s_box = ((0x0 , 0x7 , 0xA , 0xE , 0x4 , 0x1 , 0x6 , 0xD , 0x5 , 0x9 , 0x8 , 0xF , 0xC , 0x2 , 0xB, 0x3),);     # LICID	
elif (sys.argv[1] == "ULC"):
  s_box = ((0x6 , 0x5 , 0xC , 0xA , 0x1 , 0xE , 0x7 , 0x9 , 0xB , 0x0 , 0x3 , 0xD , 0x8 , 0xF , 0x4, 0x2),);     # ULC				
elif (sys.argv[1] == "SLA"):
    s_box = ((0xF , 0x8 , 0x3 , 0xE , 0x0 , 0x7 , 0xB , 0xA , 0x5 , 0xD , 0x9 , 0xC, 0x6 , 0x4, 0x2, 0x1),);     # SLA	
elif (sys.argv[1] == "PRESENT"):
  s_box = ((0xC , 0x5 , 0x6 , 0xB , 0x9 , 0x0 , 0xA , 0xD , 0x3 , 0xE , 0xF , 0x8 , 0x4 , 0x7 , 0x1, 0x2),);     # PRESENT
elif (sys.argv[1] == "WARP"):
  s_box = ((0xc,0xa, 0xd,0x3, 0xe,0xb, 0xf, 0x7, 0x8, 0x9, 0x1, 0x5, 0x0, 0x2, 0x4, 0x6),);     # WARP
elif (sys.argv[1] == "GIFT"):
  s_box = ((0x1,0xa, 0x4,0xc, 0x6,0xf, 0x3, 0x9, 0x2, 0xd, 0xb, 0x7, 0x5, 0x0, 0x8, 0xe),);       # GIFT       
elif (sys.argv[1] == "TWINE"):
  s_box = ((0xc,0x0, 0xf,0xa, 0x2,0xb, 0x9, 0x5, 0x8, 0x3, 0xd, 0x7, 0x1, 0xe, 0x6, 0x4),);       # TWINE         
elif (sys.argv[1] == "ASCON"):
  s_box = ((0x4,0xb,0x1f,0x14,0x1a,0x15,0x9,0x2,0x1b,0x5,0x8,0x12,0x1d,0x3,0x6,0x1c,), (0x1e,0x13,0x7,0xe,0x0,0xd,0x11,0x18,0x10,0xc,0x1,0x19,0x16,0xa,0xf,0x17));       # ASCON  
elif (sys.argv[1] == "FIDES-6"):
  s_box = ( (54,0,48,13,15,18,35,53,63,25,45,52,3,20,33,41),(8,10,57,37,59,36,34,2,26,50,58,24,60,19,14,42),(46,61,5,49,31,11,28,4,12,30,55,22,9,6,32,23),(27,39,21,17,16,29,62,1,40,47,51,56,7,43,38,44));#FIDES-6
elif (sys.argv[1] == "FIDES-5"):
  s_box = ((1,0,25,26,17,29,21,27,20,5,4,23,14,18,2,28),(15,8,6,3,13,7,24,16,30,9,31,10,22,12,11,19)); #FIDES-5
elif (sys.argv[1] == "SC2000-5"):
  s_box = ((20,26,7,31,19,12,10,15,22,30,13,14, 4,24, 9,18),(27,11, 1,21, 6,16, 2,28,23, 5, 8, 3, 0,17,29,25));
elif (sys.argv[1] == "SC2000-6"):
  s_box = ((47,59,25,42,15,23,28,39,26,38,36,19,60,24,29,56),(37,63,20,61,55, 2,30,44, 9,10, 6,22,53,48,51,11),(62,52,35,18,14,46, 0,54,17,40,27, 4,31, 8, 5,12),(3,16,41,34,33, 7,45,49,50,58, 1,21,43,57,32,13));
elif (sys.argv[1] == "APN-6"):
  s_box = ((0,54,48,13,15,18,53,35,25,63,45,52,3,20,41,33),(59,36,2,34,10,8,57,37,60,19,42,14,50,26,58,24),(39,27,21,17,16,29,1,62,47,40,51,56,7,43,44,38),(31,11,4,28,61,46,5,49,9,6,23,32,30,12,55,22));
elif (sys.argv[1] == "IVLBC"):
    s_box = ((0x0 , 0xF , 0xE , 0x5 , 0xD , 0x3 , 0x6 , 0xC , 0xB , 0x9 , 0xA , 0x8 , 0x7 , 0x4 , 0x2, 0x1),);
elif (sys.argv[1] == "ESLICE"):
    s_box = ((0x6 , 0x0 , 0x8 , 0xF , 0xC , 0x3 , 0x7 , 0xD , 0xB , 0xE , 0x1 , 0x4 , 0x5 , 0x9 , 0xA, 0x2),);

DDT_SIZE = (len(s_box)*len(s_box[0]))
input_size = int(math.log(DDT_SIZE,2))
DDT = np.zeros( (DDT_SIZE,DDT_SIZE) )
DDT = DDT.astype(int)
sbox_val = []


for p2 in range(DDT_SIZE):
    row = p2 >> 4
    col = p2 & 15
    sbox_val.append(s_box[row][col]);

for p1 in range(DDT_SIZE):
	for p2 in range(DDT_SIZE):
		XOR_IN = np.bitwise_xor(p1,p2);
		XOR_OUT = np.bitwise_xor(sbox_val[p1],sbox_val[p2]);
		DDT[XOR_IN][XOR_OUT] += 1


diff_arr = []
diff_arr_qm = []
diff_arr_with_1 = []
impossible_diff_arr=[]
impossible_diff_arr_qm=[]
impossible_diff_arr_new=[]
if ("sbox" in sys.argv[2]):
    check = False
    prob_point = 0
    if ("_" in sys.argv[2]):
       prob_point = int(sys.argv[2].split("_")[1])
       check = True
    for row in range(len(DDT)):
            row_hex = bin(row)[2:].zfill(input_size);
            row_arr = [int(i) for i in row_hex];
            for col in range(len(DDT[row])):
                col_hex = bin(col)[2:].zfill(input_size);
                col_arr = [int(i) for i in col_hex];
                    
                if ((DDT[row][col]==prob_point)==check):
                    diff_arr += [row_arr+col_arr];
                    diff_arr_with_1 += [[1]+row_arr+col_arr];
                    diff_arr_qm += [row_hex+col_hex];
                else:
                    impossible_diff_arr += [[1]+row_arr+col_arr];
                    impossible_diff_arr_qm += [row_hex+col_hex];
                    impossible_diff_arr_new += [row_arr+col_arr];
                

unique_entries = np.unique(DDT)
unique_entries_count = len(np.unique(DDT))
print(unique_entries)
if (sys.argv[2] == "prob"):               
    diff_arr = []
    diff_arr_with_1 = []
    impossible_diff_arr=[]
    impossible_diff_arr_qm=[]
    impossible_diff_arr_new=[]
    for row in range(len(DDT)):
            row_bin = bin(row)[2:].zfill(input_size);
            row_arr = [int(i) for i in row_bin];
            for col in range(len(DDT[row])):
                col_bin = bin(col)[2:].zfill(input_size);
                col_arr = [int(i) for i in col_bin];
                if(DDT[row][col]!=0):
                    DDT_bin = "".join(['0']*(unique_entries_count-2))
                    for num in range(1,len(unique_entries)-1):
                        if (DDT[row][col] == unique_entries[num]):
                            DDT_bin = DDT_bin[0:len(DDT_bin) - np.where(unique_entries==unique_entries[num])[0][0]] +  '1' + DDT_bin[len(DDT_bin) - np.where(unique_entries==unique_entries[num])[0][0] +1:] 
                        
                    int_DDT_bin = int(DDT_bin,2);
                    DDT_arr = [int(i) for i in DDT_bin];
                    diff_arr += [row_arr+col_arr+DDT_arr];
                    diff_arr_with_1 += [[1]+row_arr+col_arr+DDT_arr];
                    for k in range(0,pow(2,unique_entries_count-2)):
                        if (k!=int_DDT_bin):
                            im_bin = bin(k)[2:].zfill(unique_entries_count-2);
                            im_arr = [int(i) for i in im_bin];
                            impossible_diff_arr += [[1] + row_arr+col_arr+im_arr];
                            impossible_diff_arr_qm += [row_bin+col_bin+im_bin]
                else:
                    for k in range(0,pow(2,unique_entries_count-2)):
                            im_bin = bin(k)[2:].zfill(unique_entries_count-2);
                            im_arr = [int(i) for i in im_bin];
                            impossible_diff_arr += [[1]+row_arr+col_arr+im_arr];
                            impossible_diff_arr_qm += [row_bin+col_bin+im_bin];
print_DDT(DDT)
ineq_list = []
print(">>> No. of possible points: " + str(len(diff_arr_with_1)))
print(">>> No. of impossible points: " + str(len(impossible_diff_arr)))

# ########Convex Hull Linear Inequalities#################################  
vertices = map(array, diff_arr)
A, b = compute_polytope_halfspaces(vertices)

ineq_list = ineq_list + np.column_stack((b,A)).astype(int).tolist()

print("Number of Inequalities: " + str(len(ineq_list)))

if (sys.argv[4]!="-"):   
    ineq_impoints_remove = []
    for i in range(0,len(ineq_list)):
        ineq_solve_count = (np.multiply(np.array(impossible_diff_arr),np.array(ineq_list[i]))).sum(1)
        ineq_impoints_remove.append([i,len(np.where(ineq_solve_count<0)[0])])
    ineq_impoints_remove.sort(key=lambda x: x[1],reverse=True)
    final_ineq_list = np.empty((0,len(ineq_list[0])))
    offset = int(sys.argv[4])
    if (offset == 0):
        offset = len(ineq_list)
    if (int(sys.argv[5]) != 1):
        for ineq in range(0,len(ineq_list),offset):
            sub_ineq_list = [ineq_list[x[0]] for x in ineq_impoints_remove[ineq:min(ineq+offset,len(ineq_list))]]
            
            if(int(sys.argv[5])==2):
            
                for i in range(0,len(sub_ineq_list)):
                    vector1 = np.array(sub_ineq_list)
                    vector2 = np.array(sub_ineq_list[i:] + sub_ineq_list[0:i] )
                    final_ineq_list = np.vstack((final_ineq_list,(vector1 + vector2)))
            
            if(int(sys.argv[5])==3):
                for i in range(0,len(sub_ineq_list)):
                    for j in range(0,len(sub_ineq_list)):
                        vector1 = np.array(sub_ineq_list)
                        vector2 = np.array(sub_ineq_list[i:] + sub_ineq_list[0:i] )
                        vector3 = np.array(sub_ineq_list[j:] + sub_ineq_list[0:j] )
                        final_ineq_list = np.vstack((final_ineq_list,(vector1 + vector2 + vector3)))
            
            if(int(sys.argv[5])==4):
                for i in range(0,len(sub_ineq_list)):
                    for j in range(0,len(sub_ineq_list)):
                        for k in range(0,len(sub_ineq_list)):    
                            vector1 = np.array(sub_ineq_list)
                            vector2 = np.array(sub_ineq_list[i:] + sub_ineq_list[0:i] )
                            vector3 = np.array(sub_ineq_list[j:] + sub_ineq_list[0:j] )
                            vector4 = np.array(sub_ineq_list[k:] + sub_ineq_list[0:k] )
                            final_ineq_list = np.vstack((final_ineq_list,(vector1 + vector2 + vector3 + vector4)))
            if(int(sys.argv[5])==5):
                for i in range(0,len(sub_ineq_list)):
                    for j in range(0,len(sub_ineq_list)):
                        for k in range(0,len(sub_ineq_list)):    
                            for l in range(0,len(sub_ineq_list)):  
                                vector1 = np.array(sub_ineq_list)
                                vector2 = np.array(sub_ineq_list[i:] + sub_ineq_list[0:i] )
                                vector3 = np.array(sub_ineq_list[j:] + sub_ineq_list[0:j] )
                                vector4 = np.array(sub_ineq_list[k:] + sub_ineq_list[0:k] )
                                vector5 = np.array(sub_ineq_list[l:] + sub_ineq_list[0:l] )
                                final_ineq_list = np.vstack((final_ineq_list,(vector1 + vector2 + vector3 + vector4 + vector5)))
        
        ineq_list = final_ineq_list
    print(">>> Final No. of Inequalities After Addition: "+ str(len(ineq_list)))

ineq_list = MILP_Solve(np.array(ineq_list),np.array(impossible_diff_arr))

