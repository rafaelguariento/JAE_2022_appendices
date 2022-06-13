# Code for generating the results of "Prey defense phenotype mediates multiple-predator effects in tri-trophic food-webs"
# This script generate the results for beta values equal to 1

library(lattice)
library(RColorBrewer)
library(deSolve)
library(Orcs)


# Script for low risk effect

N_P1 <- 10;
N_P2 <- 10;


#vectors for storing results
result_A_r2 <- matrix(data = 0, nrow = N_P1, ncol = N_P2)
result_V1_r2 <- matrix(data = 0, nrow = N_P1, ncol = N_P2)
result_V2_r2 <- matrix(data = 0, nrow = N_P1, ncol = N_P2)
result_V12_r2 <- matrix(data = 0, nrow = N_P1, ncol = N_P2)
result_P_r2 <- matrix(data = 0, nrow = N_P1, ncol = N_P2)
result_P_no_beh_r2 <- matrix(data = 0, nrow = N_P1, ncol = N_P2)
result_density_r2 <- matrix(data = 0, nrow = N_P1, ncol = N_P2)
  

state <- c(
  A = 10,  # phenotype with no defenses
  V1 = 0,  # phenotype with defenses against predator A
  V2 = 0,  # phenotype with defenses against predator B
  V12 = 0, # phenotype with defenses against both predators
  R = 10,  # Resources under DMII + TMII
  R_asterisc = 10  # Resources under DMII
)

for (i in 1:N_P1) {
  for (j in 1:N_P2) {
    parameter <-
      c(
        f = 0.06,
        CR = 0.2,
        P1 = (i / N_P1),
        P2 = (j / N_P2),
        PR = 0.25,
        r = 1,
        KR = 100,
        C_V1 = 0.02,
        C_V2 = 0.02,
        C_V12 = 0.03,
         teta=1,
        chi_max=0.01,
        h = 0.1
      )
    
    CRM <- function(t, state, parameter) {
      
     
      
      delta_f_A_V1= parameter["chi_max"] *  (1 / (1 +exp(parameter["teta"]*((((parameter["f"])/(1+parameter["h"]*parameter["f"]) * state["R"]   )   - parameter["P1"]/(1+parameter["h"]*parameter["P1"])  - parameter["P2"]/(1+parameter["h"]*parameter["P2"]))-((((parameter["f"]-parameter["C_V1"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V1"])) * state["R"]     - parameter["CR"]))    - parameter["P2"]/(1+parameter["h"]*parameter["P2"])- parameter["PR"]*parameter["P1"]/(1+parameter["h"]*parameter["P1"]))))))
     
      delta_f_A_V2= parameter["chi_max"] *  (1 / (1 +exp(parameter["teta"]*((((parameter["f"])/(1+parameter["h"]*parameter["f"]) * state["R"]   )   - parameter["P1"]/(1+parameter["h"]*parameter["P1"])  - parameter["P2"]/(1+parameter["h"]*parameter["P2"]) )-((((parameter["f"]-parameter["C_V2"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V2"])) * state["R"]   - parameter["CR"]))   - parameter["P1"]/(1+parameter["h"]*parameter["P1"]) - parameter["PR"]*parameter["P2"]/(1+parameter["h"]*parameter["P2"]))))))
      
      delta_f_V1_V12= parameter["chi_max"] *  (1 / (1 +exp(parameter["teta"]*(((((parameter["f"]-parameter["C_V1"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V1"])) * state["R"]     - parameter["CR"]))    - parameter["P2"]/(1+parameter["h"]*parameter["P2"])- parameter["PR"]*parameter["P1"]/(1+parameter["h"]*parameter["P1"]) )-((((parameter["f"]-parameter["C_V12"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V12"])) * state["R"]    - 2 * parameter["CR"]- parameter["PR"]*parameter["P1"]/(1+parameter["h"]*parameter["P1"])- parameter["PR"]*parameter["P2"]/(1+parameter["h"]*parameter["P2"])))  ))))) 
      
      delta_f_V2_V12= parameter["chi_max"] *  (1 / (1 +exp(parameter["teta"]*(((((parameter["f"]-parameter["C_V2"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V2"])) * state["R"]   - parameter["CR"]))    - parameter["P1"]/(1+parameter["h"]*parameter["P1"])- parameter["PR"]*parameter["P2"]/(1+parameter["h"]*parameter["P2"]))-((((parameter["f"]-parameter["C_V12"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V12"])) * state["R"]    - 2 * parameter["CR"]- parameter["PR"]*parameter["P1"]/(1+parameter["h"]*parameter["P1"])- parameter["PR"]*parameter["P2"]/(1+parameter["h"]*parameter["P2"])))  ))))) 

     
      
      delta_f_V1_A= parameter["chi_max"] *  (1 / (1 +exp(parameter["teta"]*(-(((parameter["f"])/(1+parameter["h"]*parameter["f"]) * state["R"]   )   - parameter["P1"]/(1+parameter["h"]*parameter["P1"])  - parameter["P2"]/(1+parameter["h"]*parameter["P2"]))+((((parameter["f"]-parameter["C_V1"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V1"])) * state["R"]     - parameter["CR"]))    - parameter["P2"]/(1+parameter["h"]*parameter["P2"])- parameter["PR"]*parameter["P1"]/(1+parameter["h"]*parameter["P1"]))))))
      
      delta_f_V2_A= parameter["chi_max"] *  (1 / (1 +exp(parameter["teta"]*(-(((parameter["f"])/(1+parameter["h"]*parameter["f"]) * state["R"]   )   - parameter["P1"]/(1+parameter["h"]*parameter["P1"])  - parameter["P2"]/(1+parameter["h"]*parameter["P2"]) )+((((parameter["f"]/(1+parameter["h"]*(parameter["f"]-parameter["C_V2"]))-parameter["C_V2"]) * state["R"]   - parameter["CR"]))   - parameter["P1"]/(1+parameter["h"]*parameter["P1"]) - parameter["PR"]*parameter["P2"]/(1+parameter["h"]*parameter["P2"]))))))
      
      delta_f_V12_V1= parameter["chi_max"] *  (1 / (1 +exp(parameter["teta"]*(-((((parameter["f"]-parameter["C_V1"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V1"])) * state["R"]     - parameter["CR"]))    - parameter["P2"]/(1+parameter["h"]*parameter["P2"])- parameter["PR"]*parameter["P1"]/(1+parameter["h"]*parameter["P1"]) )+((((parameter["f"]-parameter["C_V12"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V12"])) * state["R"]    - 2 * parameter["CR"]- parameter["PR"]*parameter["P1"]/(1+parameter["h"]*parameter["P1"])- parameter["PR"]*parameter["P2"]/(1+parameter["h"]*parameter["P2"])))  )))))
      
      delta_f_V12_V2= parameter["chi_max"] *  (1 / (1 +exp(parameter["teta"]*(-((((parameter["f"]-parameter["C_V2"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V2"])) * state["R"]   - parameter["CR"]))    - parameter["P1"]/(1+parameter["h"]*parameter["P1"])- parameter["PR"]*parameter["P2"]/(1+parameter["h"]*parameter["P2"]))+((((parameter["f"]-parameter["C_V12"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V12"])) * state["R"]    - 2 * parameter["CR"]- parameter["PR"]*parameter["P1"]/(1+parameter["h"]*parameter["P1"])- parameter["PR"]*parameter["P2"]/(1+parameter["h"]*parameter["P2"])))  )))))
      
      
      
      
       # rates of change
      
      dA <-
        ((parameter["f"])/(1+parameter["h"]*parameter["f"]) * state["R"]   *state["A"])   - parameter["P1"]/(1+parameter["h"]*parameter["P1"]) * state["A"] - parameter["P2"]/(1+parameter["h"]*parameter["P2"]) * state["A"] - delta_f_A_V1* state["A"]-delta_f_A_V2* state["A"] + delta_f_V1_A* state["V1"] + delta_f_V2_A* state["V2"]
      dV1 <-
        (((parameter["f"]-parameter["C_V1"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V1"])) * state["R"] - parameter["CR"])*state["V1"])  + delta_f_A_V1* state["A"] - delta_f_V1_A* state["V1"] - parameter["PR"]* parameter["P1"]/(1+parameter["h"]*parameter["P1"]) * state["V1"]  - parameter["P2"]/(1+parameter["h"]*parameter["P2"]) * state["V1"] - delta_f_V1_V12* state["V1"] + delta_f_V12_V1* state["V12"]
      dV2 <-
        (((parameter["f"]-parameter["C_V2"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V2"])) * state["R"]   - parameter["CR"])*state["V2"])    + delta_f_A_V2* state["A"] - delta_f_V2_A* state["V2"] - parameter["PR"]*parameter["P2"]/(1+parameter["h"]*parameter["P2"]) * state["V2"]- parameter["P1"]/(1+parameter["h"]*parameter["P1"]) * state["V2"] - delta_f_V2_V12* state["V2"] + delta_f_V12_V2* state["V12"]
      dV12 <-
        (((parameter["f"]-parameter["C_V12"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V12"])) * state["R"]    - 2 * parameter["CR"])*state["V12"])    + delta_f_V1_V12* state["V1"] + delta_f_V2_V12* state["V2"] - delta_f_V12_V1* state["V12"] - delta_f_V12_V2* state["V12"] - parameter["PR"]*parameter["P1"]/(1+parameter["h"]*parameter["P1"]) * state["V12"]- parameter["PR"]*parameter["P2"]/(1+parameter["h"]*parameter["P2"]) * state["V12"]
      dP <-
        parameter["r"] * (1 - state["R"] / parameter["KR"]) * state["R"] - (parameter["f"])/(1+parameter["h"]*parameter["f"]) *state["R"] * state["A"]    - (parameter["f"]-parameter["C_V1"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V1"])) *state["R"] * state["V1"]     - (parameter["f"]-parameter["C_V2"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V2"])) *state["R"] * state["V2"]   - (parameter["f"]-parameter["C_V12"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V12"])) *state["R"] * state["V12"]   
      dPP <-
        parameter["r"] * (1 - state["R_asterisc"] / parameter["KR"]) * state["R_asterisc"] - (parameter["f"])/(1+parameter["h"]*parameter["f"]) *state["R_asterisc"] * state["A"]    - (parameter["f"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V1"])) *state["R_asterisc"] * state["V1"]    - (parameter["f"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V2"])) *state["R_asterisc"] * state["V2"]    - (parameter["f"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V12"])) *state["R_asterisc"] * state["V12"]   
      
      # return the rate of change
      list(c(dA, dV1, dV2, dV12, dP, dPP))
      
    }
     times <- seq(0, 2000, by = 0.5)
    
    
    result <-
      ode(
        y = state,
        times = times,
        func = CRM,
        parms = parameter
      )
    head(result)
    result_A_r2[i, j] <- result[3990, "A"]
    result_V1_r2[i, j] <- result[3990, "V1"]
    result_V2_r2[i, j] <- result[3990, "V2"]
    result_V12_r2[i, j] <- result[3990, "V12"]
    result_P_r2[i, j] <- result[3990, "R"]
    result_P_no_beh_r2[i, j] <- result[3990, "R_asterisc"]
    result_density_r2[i, j] <-
      result[3990, "A"] + result[3990, "V1"] + result[3990, "V2"] + result[3990, "V12"]
    
  }
}

# End of the script for low risk effect

# script for medium risk effect

N_P1 <- 10;
N_P2 <- 10;

result_A_r3 <- matrix(data = 0, nrow = N_P1, ncol = N_P2)
result_V1_r3 <- matrix(data = 0, nrow = N_P1, ncol = N_P2)
result_V2_r3 <- matrix(data = 0, nrow = N_P1, ncol = N_P2)
result_V12_r3 <- matrix(data = 0, nrow = N_P1, ncol = N_P2)
result_P_r3 <- matrix(data = 0, nrow = N_P1, ncol = N_P2)
result_P_no_beh_r3 <- matrix(data = 0, nrow = N_P1, ncol = N_P2)
result_density_r3 <- matrix(data = 0, nrow = N_P1, ncol = N_P2)



state <- c(
  A = 10,
  V1 = 0,
  V2 = 0,
  V12 = 0,
  R = 10,
  R_asterisc = 10
)

for (i in 1:N_P1) {
  for (j in 1:N_P2) {
    parameter <-
      c(
        f = 0.06,
         
        CR = 0.25,
        P1 = (i / N_P1),
        P2 = (j / N_P2),
        PR = 0.25,
        r = 1,
        KR = 100,
        
        C_V1 = 0.02,
        C_V2 = 0.02,
        C_V12 = 0.03,
         teta=1,
        chi_max=0.01,
        h = 0.1
      )
    
    CRM <- function(t, state, parameter) {
      
      
      
      delta_f_A_V1= parameter["chi_max"] *  (1 / (1 +exp(parameter["teta"]*((((parameter["f"])/(1+parameter["h"]*parameter["f"]) * state["R"]   )   - parameter["P1"]/(1+parameter["h"]*parameter["P1"])  - parameter["P2"]/(1+parameter["h"]*parameter["P2"]))-((((parameter["f"]-parameter["C_V1"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V1"])) * state["R"]     - parameter["CR"]))    - parameter["P2"]/(1+parameter["h"]*parameter["P2"])- parameter["PR"]*parameter["P1"]/(1+parameter["h"]*parameter["P1"]))))))
      
      delta_f_A_V2= parameter["chi_max"] *  (1 / (1 +exp(parameter["teta"]*((((parameter["f"])/(1+parameter["h"]*parameter["f"]) * state["R"]   )   - parameter["P1"]/(1+parameter["h"]*parameter["P1"])  - parameter["P2"]/(1+parameter["h"]*parameter["P2"]) )-((((parameter["f"]-parameter["C_V2"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V2"])) * state["R"]   - parameter["CR"]))   - parameter["P1"]/(1+parameter["h"]*parameter["P1"]) - parameter["PR"]*parameter["P2"]/(1+parameter["h"]*parameter["P2"]))))))
      
      delta_f_V1_V12= parameter["chi_max"] *  (1 / (1 +exp(parameter["teta"]*(((((parameter["f"]-parameter["C_V1"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V1"])) * state["R"]     - parameter["CR"]))    - parameter["P2"]/(1+parameter["h"]*parameter["P2"])- parameter["PR"]*parameter["P1"]/(1+parameter["h"]*parameter["P1"]) )-((((parameter["f"]-parameter["C_V12"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V12"])) * state["R"]    - 2 * parameter["CR"]- parameter["PR"]*parameter["P1"]/(1+parameter["h"]*parameter["P1"])- parameter["PR"]*parameter["P2"]/(1+parameter["h"]*parameter["P2"])))  ))))) 
      
      delta_f_V2_V12= parameter["chi_max"] *  (1 / (1 +exp(parameter["teta"]*(((((parameter["f"]-parameter["C_V2"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V2"])) * state["R"]   - parameter["CR"]))    - parameter["P1"]/(1+parameter["h"]*parameter["P1"])- parameter["PR"]*parameter["P2"]/(1+parameter["h"]*parameter["P2"]))-((((parameter["f"]-parameter["C_V12"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V12"])) * state["R"]    - 2 * parameter["CR"]- parameter["PR"]*parameter["P1"]/(1+parameter["h"]*parameter["P1"])- parameter["PR"]*parameter["P2"]/(1+parameter["h"]*parameter["P2"])))  ))))) 
      
      
      
      delta_f_V1_A= parameter["chi_max"] *  (1 / (1 +exp(parameter["teta"]*(-(((parameter["f"])/(1+parameter["h"]*parameter["f"]) * state["R"]   )   - parameter["P1"]/(1+parameter["h"]*parameter["P1"])  - parameter["P2"]/(1+parameter["h"]*parameter["P2"]))+((((parameter["f"]-parameter["C_V1"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V1"])) * state["R"]     - parameter["CR"]))    - parameter["P2"]/(1+parameter["h"]*parameter["P2"])- parameter["PR"]*parameter["P1"]/(1+parameter["h"]*parameter["P1"]))))))
      
      delta_f_V2_A= parameter["chi_max"] *  (1 / (1 +exp(parameter["teta"]*(-(((parameter["f"])/(1+parameter["h"]*parameter["f"]) * state["R"]   )   - parameter["P1"]/(1+parameter["h"]*parameter["P1"])  - parameter["P2"]/(1+parameter["h"]*parameter["P2"]) )+((((parameter["f"]/(1+parameter["h"]*(parameter["f"]-parameter["C_V2"]))-parameter["C_V2"]) * state["R"]   - parameter["CR"]))   - parameter["P1"]/(1+parameter["h"]*parameter["P1"]) - parameter["PR"]*parameter["P2"]/(1+parameter["h"]*parameter["P2"]))))))
      
      delta_f_V12_V1= parameter["chi_max"] *  (1 / (1 +exp(parameter["teta"]*(-((((parameter["f"]-parameter["C_V1"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V1"])) * state["R"]     - parameter["CR"]))    - parameter["P2"]/(1+parameter["h"]*parameter["P2"])- parameter["PR"]*parameter["P1"]/(1+parameter["h"]*parameter["P1"]) )+((((parameter["f"]-parameter["C_V12"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V12"])) * state["R"]    - 2 * parameter["CR"]- parameter["PR"]*parameter["P1"]/(1+parameter["h"]*parameter["P1"])- parameter["PR"]*parameter["P2"]/(1+parameter["h"]*parameter["P2"])))  )))))
      
      delta_f_V12_V2= parameter["chi_max"] *  (1 / (1 +exp(parameter["teta"]*(-((((parameter["f"]-parameter["C_V2"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V2"])) * state["R"]   - parameter["CR"]))    - parameter["P1"]/(1+parameter["h"]*parameter["P1"])- parameter["PR"]*parameter["P2"]/(1+parameter["h"]*parameter["P2"]))+((((parameter["f"]-parameter["C_V12"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V12"])) * state["R"]    - 2 * parameter["CR"]- parameter["PR"]*parameter["P1"]/(1+parameter["h"]*parameter["P1"])- parameter["PR"]*parameter["P2"]/(1+parameter["h"]*parameter["P2"])))  )))))
      
      
      
      
      # rates of change
      
      dA <-
        ((parameter["f"])/(1+parameter["h"]*parameter["f"]) * state["R"]   *state["A"])   - parameter["P1"]/(1+parameter["h"]*parameter["P1"]) * state["A"] - parameter["P2"]/(1+parameter["h"]*parameter["P2"]) * state["A"] - delta_f_A_V1* state["A"]-delta_f_A_V2* state["A"] + delta_f_V1_A* state["V1"] + delta_f_V2_A* state["V2"]
      dV1 <-
        (((parameter["f"]-parameter["C_V1"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V1"])) * state["R"] - parameter["CR"])*state["V1"])  + delta_f_A_V1* state["A"] - delta_f_V1_A* state["V1"] - parameter["PR"]* parameter["P1"]/(1+parameter["h"]*parameter["P1"]) * state["V1"]  - parameter["P2"]/(1+parameter["h"]*parameter["P2"]) * state["V1"] - delta_f_V1_V12* state["V1"] + delta_f_V12_V1* state["V12"]
      dV2 <-
        (((parameter["f"]-parameter["C_V2"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V2"])) * state["R"]   - parameter["CR"])*state["V2"])    + delta_f_A_V2* state["A"] - delta_f_V2_A* state["V2"] - parameter["PR"]*parameter["P2"]/(1+parameter["h"]*parameter["P2"]) * state["V2"]- parameter["P1"]/(1+parameter["h"]*parameter["P1"]) * state["V2"] - delta_f_V2_V12* state["V2"] + delta_f_V12_V2* state["V12"]
      dV12 <-
        (((parameter["f"]-parameter["C_V12"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V12"])) * state["R"]    - 2 * parameter["CR"])*state["V12"])    + delta_f_V1_V12* state["V1"] + delta_f_V2_V12* state["V2"] - delta_f_V12_V1* state["V12"] - delta_f_V12_V2* state["V12"] - parameter["PR"]*parameter["P1"]/(1+parameter["h"]*parameter["P1"]) * state["V12"]- parameter["PR"]*parameter["P2"]/(1+parameter["h"]*parameter["P2"]) * state["V12"]
      dP <-
        parameter["r"] * (1 - state["R"] / parameter["KR"]) * state["R"] - (parameter["f"])/(1+parameter["h"]*parameter["f"]) *state["R"] * state["A"]    - (parameter["f"]-parameter["C_V1"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V1"])) *state["R"] * state["V1"]     - (parameter["f"]-parameter["C_V2"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V2"])) *state["R"] * state["V2"]   - (parameter["f"]-parameter["C_V12"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V12"])) *state["R"] * state["V12"]   
      dPP <-
        parameter["r"] * (1 - state["R_asterisc"] / parameter["KR"]) * state["R_asterisc"] - (parameter["f"])/(1+parameter["h"]*parameter["f"]) *state["R_asterisc"] * state["A"]    - (parameter["f"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V1"])) *state["R_asterisc"] * state["V1"]    - (parameter["f"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V2"])) *state["R_asterisc"] * state["V2"]    - (parameter["f"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V12"])) *state["R_asterisc"] * state["V12"]    
      
      # return the rate of change
      list(c(dA, dV1, dV2, dV12, dP, dPP))
      
    }
     times <- seq(0, 2000, by = 0.5)
    
    
    result <-
      ode(
        y = state,
        times = times,
        func = CRM,
        parms = parameter
      )
    head(result)
    result_A_r3[i, j] <- result[3990, "A"]
    result_V1_r3[i, j] <- result[3990, "V1"]
    result_V2_r3[i, j] <- result[3990, "V2"]
    result_V12_r3[i, j] <- result[3990, "V12"]
    result_P_r3[i, j] <- result[3990, "R"]
    result_P_no_beh_r3[i, j] <- result[3990, "R_asterisc"]
    result_density_r3[i, j] <-
      result[3990, "A"] + result[3990, "V1"] + result[3990, "V2"] + result[3990, "V12"]
    
  }
}



# End of the script for medium risk effect



# script for high risk effect

N_P1 <- 10;
N_P2 <- 10;

result_A_r5 <- matrix(data = 0, nrow = N_P1, ncol = N_P2)
result_V1_r5 <- matrix(data = 0, nrow = N_P1, ncol = N_P2)
result_V2_r5 <- matrix(data = 0, nrow = N_P1, ncol = N_P2)
result_V12_r5 <- matrix(data = 0, nrow = N_P1, ncol = N_P2)
result_P_r5 <- matrix(data = 0, nrow = N_P1, ncol = N_P2)
result_P_no_beh_r5 <- matrix(data = 0, nrow = N_P1, ncol = N_P2)
result_density_r5 <- matrix(data = 0, nrow = N_P1, ncol = N_P2)



state <- c(
  A = 10,
  V1 = 0,
  V2 = 0,
  V12 = 0,
  R = 10,
  R_asterisc = 10
)

for (i in 1:N_P1) {
  for (j in 1:N_P2) {
    parameter <-
      c(
        f = 0.06,
         
        CR = 0.3,
        P1 = (i / N_P1),
        P2 = (j / N_P2),
        PR = 0.25,
        r = 1,
        KR = 100,
        
        C_V1 = 0.02,
        C_V2 = 0.02,
        C_V12 = 0.03,
         teta=1,
        chi_max=0.01,
        h = 0.1
      )
    
    CRM <- function(t, state, parameter) {
      
      
      
      delta_f_A_V1= parameter["chi_max"] *  (1 / (1 +exp(parameter["teta"]*((((parameter["f"])/(1+parameter["h"]*parameter["f"]) * state["R"]   )   - parameter["P1"]/(1+parameter["h"]*parameter["P1"])  - parameter["P2"]/(1+parameter["h"]*parameter["P2"]))-((((parameter["f"]-parameter["C_V1"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V1"])) * state["R"]     - parameter["CR"]))    - parameter["P2"]/(1+parameter["h"]*parameter["P2"])- parameter["PR"]*parameter["P1"]/(1+parameter["h"]*parameter["P1"]))))))
      
      delta_f_A_V2= parameter["chi_max"] *  (1 / (1 +exp(parameter["teta"]*((((parameter["f"])/(1+parameter["h"]*parameter["f"]) * state["R"]   )   - parameter["P1"]/(1+parameter["h"]*parameter["P1"])  - parameter["P2"]/(1+parameter["h"]*parameter["P2"]) )-((((parameter["f"]-parameter["C_V2"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V2"])) * state["R"]   - parameter["CR"]))   - parameter["P1"]/(1+parameter["h"]*parameter["P1"]) - parameter["PR"]*parameter["P2"]/(1+parameter["h"]*parameter["P2"]))))))
      
      delta_f_V1_V12= parameter["chi_max"] *  (1 / (1 +exp(parameter["teta"]*(((((parameter["f"]-parameter["C_V1"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V1"])) * state["R"]     - parameter["CR"]))    - parameter["P2"]/(1+parameter["h"]*parameter["P2"])- parameter["PR"]*parameter["P1"]/(1+parameter["h"]*parameter["P1"]) )-((((parameter["f"]-parameter["C_V12"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V12"])) * state["R"]    - 2 * parameter["CR"]- parameter["PR"]*parameter["P1"]/(1+parameter["h"]*parameter["P1"])- parameter["PR"]*parameter["P2"]/(1+parameter["h"]*parameter["P2"])))  ))))) 
      
      delta_f_V2_V12= parameter["chi_max"] *  (1 / (1 +exp(parameter["teta"]*(((((parameter["f"]-parameter["C_V2"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V2"])) * state["R"]   - parameter["CR"]))    - parameter["P1"]/(1+parameter["h"]*parameter["P1"])- parameter["PR"]*parameter["P2"]/(1+parameter["h"]*parameter["P2"]))-((((parameter["f"]-parameter["C_V12"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V12"])) * state["R"]    - 2 * parameter["CR"]- parameter["PR"]*parameter["P1"]/(1+parameter["h"]*parameter["P1"])- parameter["PR"]*parameter["P2"]/(1+parameter["h"]*parameter["P2"])))  ))))) 
      
      
      
      delta_f_V1_A= parameter["chi_max"] *  (1 / (1 +exp(parameter["teta"]*(-(((parameter["f"])/(1+parameter["h"]*parameter["f"]) * state["R"]   )   - parameter["P1"]/(1+parameter["h"]*parameter["P1"])  - parameter["P2"]/(1+parameter["h"]*parameter["P2"]))+((((parameter["f"]-parameter["C_V1"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V1"])) * state["R"]     - parameter["CR"]))    - parameter["P2"]/(1+parameter["h"]*parameter["P2"])- parameter["PR"]*parameter["P1"]/(1+parameter["h"]*parameter["P1"]))))))
      
      delta_f_V2_A= parameter["chi_max"] *  (1 / (1 +exp(parameter["teta"]*(-(((parameter["f"])/(1+parameter["h"]*parameter["f"]) * state["R"]   )   - parameter["P1"]/(1+parameter["h"]*parameter["P1"])  - parameter["P2"]/(1+parameter["h"]*parameter["P2"]) )+((((parameter["f"]/(1+parameter["h"]*(parameter["f"]-parameter["C_V2"]))-parameter["C_V2"]) * state["R"]   - parameter["CR"]))   - parameter["P1"]/(1+parameter["h"]*parameter["P1"]) - parameter["PR"]*parameter["P2"]/(1+parameter["h"]*parameter["P2"]))))))
      
      delta_f_V12_V1= parameter["chi_max"] *  (1 / (1 +exp(parameter["teta"]*(-((((parameter["f"]-parameter["C_V1"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V1"])) * state["R"]     - parameter["CR"]))    - parameter["P2"]/(1+parameter["h"]*parameter["P2"])- parameter["PR"]*parameter["P1"]/(1+parameter["h"]*parameter["P1"]) )+((((parameter["f"]-parameter["C_V12"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V12"])) * state["R"]    - 2 * parameter["CR"]- parameter["PR"]*parameter["P1"]/(1+parameter["h"]*parameter["P1"])- parameter["PR"]*parameter["P2"]/(1+parameter["h"]*parameter["P2"])))  )))))
      
      delta_f_V12_V2= parameter["chi_max"] *  (1 / (1 +exp(parameter["teta"]*(-((((parameter["f"]-parameter["C_V2"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V2"])) * state["R"]   - parameter["CR"]))    - parameter["P1"]/(1+parameter["h"]*parameter["P1"])- parameter["PR"]*parameter["P2"]/(1+parameter["h"]*parameter["P2"]))+((((parameter["f"]-parameter["C_V12"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V12"])) * state["R"]    - 2 * parameter["CR"]- parameter["PR"]*parameter["P1"]/(1+parameter["h"]*parameter["P1"])- parameter["PR"]*parameter["P2"]/(1+parameter["h"]*parameter["P2"])))  )))))
      
      
      
      
      # rates of change
      
      dA <-
        ((parameter["f"])/(1+parameter["h"]*parameter["f"]) * state["R"]   *state["A"])   - parameter["P1"]/(1+parameter["h"]*parameter["P1"]) * state["A"] - parameter["P2"]/(1+parameter["h"]*parameter["P2"]) * state["A"] - delta_f_A_V1* state["A"]-delta_f_A_V2* state["A"] + delta_f_V1_A* state["V1"] + delta_f_V2_A* state["V2"]
      dV1 <-
        (((parameter["f"]-parameter["C_V1"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V1"])) * state["R"] - parameter["CR"])*state["V1"])  + delta_f_A_V1* state["A"] - delta_f_V1_A* state["V1"] - parameter["PR"]* parameter["P1"]/(1+parameter["h"]*parameter["P1"]) * state["V1"]  - parameter["P2"]/(1+parameter["h"]*parameter["P2"]) * state["V1"] - delta_f_V1_V12* state["V1"] + delta_f_V12_V1* state["V12"]
      dV2 <-
        (((parameter["f"]-parameter["C_V2"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V2"])) * state["R"]   - parameter["CR"])*state["V2"])    + delta_f_A_V2* state["A"] - delta_f_V2_A* state["V2"] - parameter["PR"]*parameter["P2"]/(1+parameter["h"]*parameter["P2"]) * state["V2"]- parameter["P1"]/(1+parameter["h"]*parameter["P1"]) * state["V2"] - delta_f_V2_V12* state["V2"] + delta_f_V12_V2* state["V12"]
      dV12 <-
        (((parameter["f"]-parameter["C_V12"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V12"])) * state["R"]    - 2 * parameter["CR"])*state["V12"])    + delta_f_V1_V12* state["V1"] + delta_f_V2_V12* state["V2"] - delta_f_V12_V1* state["V12"] - delta_f_V12_V2* state["V12"] - parameter["PR"]*parameter["P1"]/(1+parameter["h"]*parameter["P1"]) * state["V12"]- parameter["PR"]*parameter["P2"]/(1+parameter["h"]*parameter["P2"]) * state["V12"]
      dP <-
        parameter["r"] * (1 - state["R"] / parameter["KR"]) * state["R"] - (parameter["f"])/(1+parameter["h"]*parameter["f"]) *state["R"] * state["A"]    - (parameter["f"]-parameter["C_V1"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V1"])) *state["R"] * state["V1"]     - (parameter["f"]-parameter["C_V2"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V2"])) *state["R"] * state["V2"]   - (parameter["f"]-parameter["C_V12"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V12"])) *state["R"] * state["V12"]   
      dPP <-
        parameter["r"] * (1 - state["R_asterisc"] / parameter["KR"]) * state["R_asterisc"] - (parameter["f"])/(1+parameter["h"]*parameter["f"]) *state["R_asterisc"] * state["A"]    - (parameter["f"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V1"])) *state["R_asterisc"] * state["V1"]    - (parameter["f"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V2"])) *state["R_asterisc"] * state["V2"]    - (parameter["f"])/(1+parameter["h"]*(parameter["f"]-parameter["C_V12"])) *state["R_asterisc"] * state["V12"]    
      
      # return the rate of change
      list(c(dA, dV1, dV2, dV12, dP, dPP))
      
    }
     times <- seq(0, 2000, by = 0.5)
    
    
    result <-
      ode(
        y = state,
        times = times,
        func = CRM,
        parms = parameter
      )
    head(result)
    result_A_r5[i, j] <- result[3990, "A"]
    result_V1_r5[i, j] <- result[3990, "V1"]
    result_V2_r5[i, j] <- result[3990, "V2"]
    result_V12_r5[i, j] <- result[3990, "V12"]
    result_P_r5[i, j] <- result[3990, "R"]
    result_P_no_beh_r5[i, j] <- result[3990, "R_asterisc"]
    result_density_r5[i, j] <-
      result[3990, "A"] + result[3990, "V1"] + result[3990, "V2"] + result[3990, "V12"]
    
  }
}

# End of the script for high risk effect

# graphs

cols <-
  (colorRampPalette(brewer.pal(6, "YlOrRd"))) 


#graphs for phenotypes

# Find the min and max values

my.min <- min(0)
my.max <- max(30)

my.at <- seq(my.min-1, my.max + 1)

p1_2_p <- levelplot(result_A_r2, at = my.at, col.regions = cols)
p2_2_p <- levelplot(result_V1_r2, at = my.at, col.regions = cols)
p3_2_p <- levelplot(result_V2_r2, at = my.at, col.regions = cols)
p4_2_p <- levelplot(result_V12_r2, at = my.at, col.regions = cols)
p1_3_p <- levelplot(result_A_r3, at = my.at, col.regions = cols)
p2_3_p <- levelplot(result_V1_r3, at = my.at, col.regions = cols)
p3_3_p <- levelplot(result_V2_r3, at = my.at, col.regions = cols)
p4_3_p <- levelplot(result_V12_r3, at = my.at, col.regions = cols)
p1_5_p <- levelplot(result_A_r5, at = my.at, col.regions = cols)
p2_5_p <- levelplot(result_V1_r5, at = my.at, col.regions = cols)
p3_5_p <- levelplot(result_V2_r5, at = my.at, col.regions = cols)
p4_5_p <- levelplot(result_V12_r5, at = my.at, col.regions = cols)

dat <- list(p1_2_p, p2_2_p, p3_2_p, p4_2_p,p1_3_p, p2_3_p, p3_3_p, p4_3_p,p1_5_p, p2_5_p, p3_5_p, p4_5_p)
p <-
  latticeCombineGrid(
    dat,  between = list(y = 1, x = 1),
    layout = c(4, 3),
    main = "",
    xlab = "Predator B Attack Rate",
    ylab = "Predator A Attack Rate"
  )

# Please, mind that in the graph that will be generated below rows depict risk effect levels (low, intermediate and high for top to bottom) and columns depict prey with no defenses, defenses against predator A, defenses against predator B and defenses against both predators, respectively.  
print(p)


#graphs for density

# Find the min and max values
my.min <- min(0)
my.max <- max(20)

my.at <- seq(my.min-1, my.max + 1)

p_2_d <- levelplot(result_density_r2, at = my.at, col.regions = cols)
p_3_d <- levelplot(result_density_r3, at = my.at, col.regions = cols)
p_5_d <- levelplot(result_density_r5, at = my.at, col.regions = cols)

dat <- list(p_2_d, p_3_d,p_5_d)
p <-
  latticeCombineGrid(
    dat,between = list(y = 1, x = 1),
    layout = c(1, 3),
    main = "",
    xlab = "Predator B Attack Rate",
    ylab = "Predator A Attack Rate"
  )

# Please, mind that in the graph that will be generated below rows depict risk effect levels (low, intermediate and high for top to bottom)
print(p)


#graphs for resources


#graphs for resources with DMII + TMII


my.min <- min(0)
my.max <- max(40)

my.at <- seq(my.min-1, my.max + 1)

p1_2_r <- levelplot(result_P_r2, at = my.at, col.regions = cols)
p1_3_r <- levelplot(result_P_r3, at = my.at, col.regions = cols)
p1_5_r <- levelplot(result_P_r5, at = my.at, col.regions = cols)

dat <- list(p1_2_r,p1_3_r,p1_5_r)
p <-
  latticeCombineGrid(
    dat,between = list(y = 1, x = 1),
    layout = c(1, 3),
    main = "",
    xlab = "Predator B Attack Rate",
    ylab = "Predator A Attack Rate"
  )

# Please, mind that in the graph that will be generated below rows depict risk effect levels (low, intermediate and high for top to bottom)
print(p)


#graphs for resources with DMII 


my.min <- min(0)
my.max <- max(40)

my.at <- seq(my.min-1, my.max + 1)


p2_2_r <- levelplot(result_P_no_beh_r2, at = my.at, col.regions = cols)

p2_3_r <- levelplot(result_P_no_beh_r3, at = my.at, col.regions = cols)

p2_5_r <- levelplot(result_P_no_beh_r5, at = my.at, col.regions = cols)

dat <- list(p2_2_r, p2_3_r, p2_5_r)
p <-
  latticeCombineGrid(
    dat,between = list(y = 1, x = 1),
    layout = c(1, 3),
    main = "",
    xlab = "Predator B Attack Rate",
    ylab = "Predator A Attack Rate"
  )

# Please, mind that in the graph that will be generated below rows depict risk effect levels (low, intermediate and high for top to bottom)   
print(p)
