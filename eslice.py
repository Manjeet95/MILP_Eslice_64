'''
Arguments:
1. Block Size - 64/128 depdends on ESLICE block size    
2. Number for rounds
3. Minimum Number of active S-boxes
4. Maximum Number of active S-boxes in each round
5. No. of trails to find
6. Solver to be used (GUROBI/CPLEX)
(if you are changing the code for another cipher then please change no. of ineq. in line 247 and 382)
python eslice.py 64 1 1 2 1 GUROBI 
python eslice.py 64 5 1 2 1 CPLEX
'''
##########################
from numpy import *
import numpy as np
from gurobipy import *
from itertools import combinations
import time
import copy
import math
import os
import shutil
#########################


import string
import sys
from math import floor

GUROBI_EXISTS = False
CPLEX_EXISTS = False

try:
  import google.colab
  IN_COLAB = True
except:
  IN_COLAB = False

if (sys.argv[6] == "GUROBI"):
    try: 
    	from gurobipy import *
    	GUROBI_EXISTS = True
    except:
    	GUROBI_EXISTS = False

if (sys.argv[6] == "CPLEX"):
	try: 
		from docplex.mp.model_reader import ModelReader
		CPLEX_EXISTS = True
	except:
		CPLEX_EXISTS = False

if (GUROBI_EXISTS == False) and (CPLEX_EXISTS == False):
	sys.exit()

P_left = (13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28,
          1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 29, 30, 31, 32)
ROUND = int(sys.argv[2])
min_sbox = int(sys.argv[3])
max_sbox_round = int(sys.argv[4])

block_size = 64
sbox_size = 4
nibble = 8
sbox = [[6, 0, 8, 15, 12, 3, 7, 13, 11, 14, 1, 4, 5, 9, 10, 2]]
var_and_num_AS = {"x":[block_size,block_size],"A":[nibble,0],"d":[block_size/2,0]}
var_and_num_DC = {"x":[block_size,block_size],"p":[nibble*2,0],"d":[block_size/2,0]}




ESLICE = int(sys.argv[1])
num_of_sbox = int(ESLICE/8)
if (ESLICE!=64):
    print("Incorrect Parameters!")
    sys.exit()
    


def obj_fun(model_goal, model_filename, BanList):
    with open(model_filename, "a") as f:
        if model_goal == "AS":
            for i in range(1, ROUND * nibble):
                f.write("A%d + " % i)
            f.write("A%d" % (ROUND * nibble))
        elif model_goal == "DC":
            for i in range(1, ROUND * nibble):
                if [i] in BanList:
                    f.write("+ 2 p%d " % (2*i-1))
                    f.write("+ 3 p%d " % (2*i))
            if [ROUND * nibble] in BanList:
                f.write("+ 2 p%d " % (2*ROUND*nibble-1))
                f.write("+ 3 p%d" % (2*ROUND*nibble))
                                
                                
def gen_input_state():
    state = np.arange(1,block_size+1)
    return state

def get_state_through_sbox(r):
    state = np.arange(block_size*r+1,block_size*r+1+block_size/2)
    return state

def get_state_through_per_right(state):
    state1 = zeros(int(block_size/2), dtype = int16)
    for i in range(int(block_size/2)):
        state1[i] = state[(i+20)%32]
    return state1

def get_state_through_per_left(state):
    state1 = np.arange(1,block_size/2+1)
    for i in range(int(block_size/2)):
        state1[i] = state[P_left[i]-1]
    return state1

def diff_propagation_of_sbox_AS(sbox_size,model_filename,var,ine):
	vector_size = 2 * sbox_size 
	symbal = [" " for i in range(vector_size+1)]
	for i in range(sbox_size-1):
		with open(model_filename, "a") as f:
			f.write("x%d + "%(var["x"][i]))
	with open(model_filename, "a") as f:
		f.write("x%d - A%d >= 0\n"%(var["x"][sbox_size-1],var["A"]))
	for i in range(sbox_size):
		with open(model_filename, "a") as f:
			f.write("A%d - x%d >= 0\n"%(var["A"],var["x"][i]))
	for i in range(len(ine)):
		for k in range (0,vector_size):
			if ine[i,k] < 0:
				symbal[k] = '-'
			elif ine[i,k] >= 0:
				symbal[k] = '+'
		for k1 in range(sbox_size):
			with open(model_filename, "a") as f:
				f.write("%s %d x%d "%(symbal[k1], abs(ine[i,k1]),var["x"][k1]))
		for k2 in range (sbox_size,sbox_size*2):
			with open(model_filename, "a") as f:
				f.write("%s %d x%d "%(symbal[k2], abs(ine[i,k2]),var["y"][k2-sbox_size]))
		if ine[i,vector_size] > 0:
			symbal[vector_size] = '-'
		elif ine[i,vector_size] <= 0:
			symbal[vector_size] = '+'
		with open(model_filename, "a") as f:
			f.write(" >= %s %d \n"%(symbal[vector_size], abs(ine[i,vector_size])))
                        
def diff_propagation_of_sbox_DC(sbox_size,num_of_p_var,model_filename,var,ine):
	vector_size = 2 * sbox_size + num_of_p_var 
	symbal = [" " for i in range(int(vector_size)+1)]
	for i in range(len(ine)):
		for k in range (0,int(vector_size)):
			if ine[i,k] < 0:
				symbal[k] = '-'
			elif ine[i,k] >= 0:
				symbal[k] = '+'
		for k1 in range(sbox_size):
			with open(model_filename, "a") as f:
				f.write("%s %d x%d "%(symbal[k1], abs(ine[i,k1]),var["x"][k1]))
		for k2 in range (sbox_size,sbox_size*2):
			with open(model_filename, "a") as f:
				f.write("%s %d x%d "%(symbal[k2], abs(ine[i,k2]),var["y"][k2-sbox_size]))
		for k3 in range(sbox_size*2,int(vector_size)):
			with open(model_filename, "a") as f:
				f.write("%s %d p%d "%(symbal[k3], abs(ine[i,k3]),var["p"][k3-sbox_size*2]))
		if ine[i,int(vector_size)] > 0:
			symbal[vector_size] = '-'
		elif ine[i,int(vector_size)] <= 0:
			symbal[int(vector_size)] = '+'
		with open(model_filename, "a") as f:
			f.write(" >= %s %d \n"%(symbal[int(vector_size)], abs(ine[i,int(vector_size)])))

def constraints_sbox(model_filename,var):
    with open(model_filename, "a") as f:
        for i in range(sbox_size):
            f.write(" 4 x%d "%(var['x'][i]))
            if i != sbox_size-1:
                f.write("+")
        for i in range(sbox_size):
            f.write(" - x%d "%(var['y'][i]))
        f.write(" >= 0\n")
        for i in range(sbox_size):
            f.write(" 4 x%d "%(var['y'][i]))
            if i != sbox_size-1:
                f.write("+")
        for i in range(sbox_size):
            f.write(" - x%d "%(var['x'][i]))
        f.write(" >= 0\n")
        
def diff_propagation_of_sbox(model_goal, model_filename, r, state, state_through_sbox, BanList):
    var = {}
    for n in range(nibble):
        var["x"] = state[sbox_size*n : sbox_size*(n+1)]
        var["y"] = state_through_sbox[sbox_size*n : sbox_size*(n+1)]
        
        if model_goal == "AS":
            var["A"] = nibble*(r-1) + n + 1
            ine = np.loadtxt("txt_AS_sbox_reduce_ineq.txt", int)
            diff_propagation_of_sbox_AS(sbox_size, model_filename, var, ine)
        
        elif model_goal == "DC":
            if [nibble*(r-1)+n+1] in BanList:
                constraints_sbox(model_filename,var)
                var["p"] = [2*n + (r-1)*2*nibble + 1, 2*n + (r-1)*2*nibble + 2]
                ine = np.loadtxt("txt_DC_sbox_reduce_ineq.txt", int)
                diff_propagation_of_sbox_DC(sbox_size, var_and_num_DC["p"][0]//nibble, model_filename, var, ine)
            elif [nibble*(r-1)+n+1] not in BanList:
                with open(model_filename, 'a') as f:
                    for i in range(sbox_size):
                        f.write("x%d" % (var["x"][i]))
                        f.write(" = 0\n")
                    for i in range(sbox_size):
                        f.write("x%d" % (var["y"][i]))
                        f.write(" = 0\n")

def get_state_through_xor(r):
    state1 = np.arange(r*block_size+block_size/2+1,r*block_size+block_size/2+1+block_size/2)
    return state1

def diff_propagation_of_xor(model_filename,r,state1,state2,state3):
    var={}
    for i in range(int(block_size/2)):
        var["x"] = state1[i]
        var["y"] = state2[i]
        var["z"] = state3[i]
        var["d"] = block_size/2*(r-1)+i+1
        diff_propagation_of_xor_bit(model_filename,var)

def diff_propagation_of_xor_bit(model_filename,var):
    with open(model_filename, "a") as f:
        f.write("x%d + x%d + x%d - 2 d%d >= 0\n"%(var["x"], var["y"], var["z"], var["d"]))
        f.write("d%d - x%d >= 0\n"%(var["d"],var["x"]))
        f.write("d%d - x%d >= 0\n"%(var["d"],var["y"]))
        f.write("d%d - x%d >= 0\n"%(var["d"],var["z"]))
        f.write("x%d + x%d + x%d <= 2\n"%(var["x"], var["y"], var["z"]))

		


def input_non_zero(model_filename):
    for i in range (1,block_size):
        with open(model_filename, "a") as f:
            f.write("x%d + "%(i))
    with open(model_filename, "a") as f:
        f.write("x%d >= 1\n"%(block_size))

def constraints_max_sboxes_round(model_filename):
    with open(model_filename, "a") as f:
        for i in range(1, ROUND + 1):
            for j in range(1, num_of_sbox + 1):
                f.write("A" + str(8 * (i - 1) + j))
                if j != num_of_sbox:
                    f.write(" + ")
            f.write(" <= " + str(max_sbox_round) + "\n")

def constraints_min_sboxes(model_filename):
    with open(model_filename, "a") as f:
        for i in range(1, ROUND + 1):
            for j in range(1, num_of_sbox + 1):
                f.write("A" + str(8 * (i - 1) + j))
                if i != ROUND or j != num_of_sbox:
                    f.write(" + ")
                else:
                    f.write(" >= " + str(min_sbox) + "\n")

def binary(model_goal, model_filename, BanList): 
    with open(model_filename, "a") as f:
        f.write("Binary\n")
        
    if model_goal == "AS":	
        var_dict = var_and_num_AS.copy()
        var_key = var_and_num_AS.keys()
    elif model_goal == "DC":	
        var_dict = var_and_num_DC.copy()
        var_key = var_and_num_DC.keys()

    for i in var_key:
        if i != "p":
            for j in range(1, int(ROUND * var_dict[i][0]) + int(var_dict[i][1]) + 1):
                with open(model_filename, "a") as f:
                    f.write(i + "%d\n" % j)
    
    for k in BanList:
        with open(model_filename, "a") as f:
            f.write("p%d\n" % (2 * k[0] - 1))
            f.write("p%d\n" % (2 * k[0]))
    
    with open(model_filename, "a") as f:
        f.write("End")


def constraints_state_feistel(model_goal,model_filename,BanList):
    state = gen_input_state()
    state_l = state[0:len(state)//2]
    state_r = state[len(state)//2:len(state)]
    for r in range (1,ROUND+1):
        sc_state_l = get_state_through_sbox(r)
        ls_state_r = get_state_through_per_right(state_r)
        p_state_l = get_state_through_per_left(sc_state_l)
        diff_propagation_of_sbox(model_goal,model_filename,r,state_l,sc_state_l,BanList)
        if r < ROUND:
            state_r = state_l
            xor_state = get_state_through_xor(r)
            diff_propagation_of_xor(model_filename,r,p_state_l,ls_state_r,xor_state)
            state_l = xor_state
        elif r == ROUND:
            xor_state = get_state_through_xor(r)
            diff_propagation_of_xor(model_filename,r,p_state_l,ls_state_r,xor_state)
            state_r = xor_state

def input_sbox_non_zero(model_filename,FixList,BanList):
    if len(FixList) == 0:
        with open(model_filename, "a") as f:
            for i in range(1,nibble+1):
                if [i] in BanList:
                    for j in range(1,sbox_size+1):
                        f.write("x" + str(4*(i-1)+j))
                        if j != sbox_size:
                            f.write(" + ")
                    f.write(" >= 1\n")


    
	    
def PrintOuter(FixList,BanList):
    model_filename = "Outer" +"_ESLICE_" + str(ESLICE) + "_" + str(ROUND) +".lp"
    model_goal="AS"
    with open(model_filename, "w") as f:
        f.write("Minimize\n")
    obj_fun(model_goal,model_filename,BanList)
    with open(model_filename, "a") as f:
        f.write("\nSubject to\n")
    constraints_state_feistel(model_goal,model_filename,BanList)
    constraints_max_sboxes_round(model_filename)
    constraints_min_sboxes(model_filename)
    input_non_zero(model_filename)
    binary(model_goal,model_filename,BanList)

def PrintInner(FixList,BanList):
    model_goal = "DC"
    model_filename = "Inner" +"_ESLICE_" + str(ESLICE) + "_" + str(ROUND) + ".lp"
    with open(model_filename, "w") as f:
        f.write("Minimize\n")
    obj_fun(model_goal,model_filename,BanList)
    with open(model_filename, "a") as f:
        f.write("\nSubject to\n")                             
    constraints_state_feistel(model_goal,model_filename,BanList)
    input_sbox_non_zero(model_filename,FixList,BanList)
    binary(model_goal,model_filename,BanList)
    #################

def read_txt(filename):
	line = []
	with open(filename) as f:
		for j in f.readlines():
			line.append(j)
	return line

def get_trail_feistel(model_goal,ROUND,count):
	if model_goal == "AS":
		x_var = [0 for j in range (ROUND*var_and_num_AS["x"][0]+var_and_num_AS["x"][1])]
	elif model_goal == "DC":
		x_var = [0 for j in range (int(ROUND*var_and_num_DC["x"][0]+var_and_num_DC["x"][1]))]	
	shutil.copy("ESLICE_" + str(ROUND) + "_"+ model_goal + "_" + str(count) + ".txt", "result_ESLICE_"+model_goal+"_"+str(ROUND)+".txt")
	var = read_txt("ESLICE_" + str(ROUND) + "_" + model_goal + "_" + str(count) + ".txt")
	
	for i in var:
		if i[0] == "x":
			x_index = int(re.findall(r'(-?[\d]+)',i)[0])
			x_value = int(re.findall(r'(-?[\d]+)',i)[1])
			x_var[x_index-1] = x_value
	file_path = "result_ESLICE_"+model_goal+"_"+str(ROUND)+"_min_Sbox_"+str(min_sbox)+"_max_S_boxes_"+str(max_sbox_round)+"_round_best_path.txt"
	with open(file_path,"w") as f:
		f.write("x_var list is " + str(x_var)+"\n")
	state = gen_input_state()
	state_l = state[0:int(len(state)/2)]
	state_r = state[int(len(state)/2):len(state)]
	for r in range(1,ROUND+1):
		with open(file_path,"a") as f:
			f.write("\n%dth round left: input diff of Sbox: \n"%(r))
		for i in range(len(state_l)):
			if i%4 == 0:
				with open(file_path,"a") as f:
					f.write(" ")
			with open(file_path,"a") as f:
				f.write(str(x_var[int(state_l[i])-1]))
		sc_state_l = get_state_through_sbox(r)
		with open(file_path,"a") as f:
			f.write("\n%dth round left: output diff of Sbox: \n"%(r))
		for i in range(len(sc_state_l)):
			if i%4 == 0:
				with open(file_path,"a") as f:
					f.write(" ")
			with open(file_path,"a") as f:
				f.write(str(x_var[int(sc_state_l[i])-1]))
		with open(file_path,"a") as f:
			f.write("\n%dth round right: \n"%(r))
		for i in range(len(state_r)):
			if i%4 == 0:
				with open(file_path,"a") as f:
					f.write(" ")
			with open(file_path,"a") as f:
				f.write(str(x_var[int(state_r[i])-1]))
		
		state_r = state_l
		xor_state = get_state_through_xor(r)
		state_l = xor_state
                        
def strtoint(s):
    s1 = ''
    result = []
    for i in range(0,len(s)):
        if s[i] >= '0' and s[i]<= '9':
            s1 = s1 + s[i]
    result.append(int(s1))  
    return result

def solve_model():
    count = 1
    fsl = []
    fslstring = []
    BanList = []
    bl = []
    FixList = []
	
    while (count <= int(sys.argv[5])):
        PrintOuter(FixList, BanList)
	    
        time_start = time.time()	
        
        if (GUROBI_EXISTS == True):
            o = read("Outer" + "_ESLICE_" + str(ESLICE) + "_" + str(ROUND) + ".lp")
            o.optimize()
            obj = o.getObjective()
        
        ##########################CPLEX####################################
        if (CPLEX_EXISTS == True):
            mr = ModelReader()
            o = mr.read("Outer" + "_ESLICE_" + str(ESLICE) + "_" + str(ROUND) + ".lp")
            o_sol = o.solve(log_output=True)
            
        ###########################GUROBI#################################
        if (GUROBI_EXISTS == True):
            o_obj = obj.getValue()
            
        ##################CPELX###############
        if (CPLEX_EXISTS == True):
            o_obj = o_sol.get_objective_value()
	    
        time_end = time.time()
        timespend = time_end - time_start
		
        optimal_solution_AS = "ESLICE_" + str(ROUND) + "_AS_" + str(count) + ".txt"									
        with open(optimal_solution_AS, "w") as f:
            f.write("solving the model " + "Outer" + "_ESLICE_" + str(ESLICE) + "_" + str(ROUND) + ".lp" + "\n")
            f.write('obj is: %d\n' % (o_obj))
            f.write('time is %d s.\n' % (timespend))
            for v in o.getVars():
                f.write("%s = %d\n" % (v.varName, int(round(v.x))))
        


        if o_obj < (min_sbox + 64) :
            b1=[]
            fsl = []
            fslstring = []
            
            ###########################GUROBI#############################
            if(GUROBI_EXISTS == True):
                for v in o.getVars():
                    if v.x == 1 and v.VarName[0] == 'A':
                        fslstring.append(v.VarName)
            ###########CPLEX###########################
            if(CPLEX_EXISTS == True):
                for A in o_sol.iter_variables():
                    if ('A' in str(v)):
                        fslstring.append(str(v))
            for f in fslstring:
                fsl.append(strtoint(f))
            for f in fslstring:
                bl.append(strtoint(f))
            BanList.append(bl)
            print("*\n*\n*\n*\n")
            print(BanList)
            print("*\n*\n*\n*\n")
            
            print(fsl)
            PrintInner(FixList,fsl)

            FixList = []
            time_start = time.time()
            ###########################GUROBI#################################
            if(GUROBI_EXISTS == True):
                i = read("Inner" +"_ESLICE_" + str(ESLICE) + "_" + str(ROUND) + ".lp")
                i.optimize()
                i_obj = i.getObjective().getValue()
            ###################CPLEX#############################
            if(CPLEX_EXISTS == True):
                i = mr.read("Inner" +"_ESLICE_" + str(ESLICE) + "_" + str(ROUND) + ".lp")
                i_sol  = i.solve(log_output=True)
                i_obj = i_sol.get_objective_value()
		
            time_end = time.time()
            timespend = time_end - time_start

            optimal_solution_AS = "ESLICE_" + str(ROUND) + "_DC_" + str(count) + ".txt"									
            with open(optimal_solution_AS, "w") as f:
                f.write("solving the model " + "Inner" + "_ESLICE_" + str(ESLICE) + "_" + str(ROUND) + ".lp" + "\n")
                f.write('obj is: %d\n' % (i_obj))
                f.write('time is %d s.\n' % (timespend))
                for v in i.getVars():
                    f.write("%s = %d\n" % (v.varName, int(round(v.x))))

            print("Number of Active S-boxes: " + str(o_obj))  
            print("FOUND Optimal Probability: " + str(i_obj))

            optimal_solution_sol="ESLICE_" + str(ROUND) + "_optimal_solution_" + str(count) + ".txt"
            with open(optimal_solution_sol, "w") as f:
                f.write(str(fsl) + "minimum weight = " + str(i_obj) + "\n")
            # get_trail_feistel('AS',ROUND,count)
            get_trail_feistel('DC',ROUND,count)
        count = count + 1

solve_model()


