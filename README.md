# MILP_Eslice_64

MILP-Based Differential Characteristics Search Problem for Block Cipher Eslice_64

## Source Code

There are two files in this source code:
- Ineq_Reduction.py
- eslice.py

## Generation and Reduction of Linear Inequalities

- First, Ineq_Reduction.py generates linear inequalities. Then, minimize the linear inequalities using GUROBI (to reduce the number of active S-boxes or to optimize the probability of differential characteristics).

- The command 'python Ineq_Reduction.py ESLICE sbox GUROBI -' gives the list of minimized linear inequalities to reduce the number of active S-boxes.

- The command 'python Ineq_Reduction.py ESLICE prob GUROBI -' provides the list of minimized linear inequalities to optimize the probability of differential characteristics.

## MILP model to minimize the number of active S-boxes and optimize the probability of differential characteristics

- Number of active S-boxes and high probability differential characteristic for 15-round Eslice_64 is searched using the following command:
```bash
    python eslice.py 64 15 5 4 1 GUROBI
```
