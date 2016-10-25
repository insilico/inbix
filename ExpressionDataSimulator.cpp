/* 
 * File:   ExpressionDataSimulator.cpp
 * Author: bwhite
 * 
 * Created on October 18, 2016, 8:02 PM
 */

#include <string>
#include <vector>
#include <random>

#include <armadillo>

#include "plink.h"
#include "options.h"

#include "ExpressionDataSimulator.h"

using namespace std;
using namespace arma;

ExpressionDataSimulator::ExpressionDataSimulator(Plink* plinkPtr):
  PP(plinkPtr) {
}

ExpressionDataSimulator::ExpressionDataSimulator(const ExpressionDataSimulator& orig) {
}

ExpressionDataSimulator::~ExpressionDataSimulator() {
}

// ----------------------------------------------------------------------------
// Create a differentially coexpressed data set without main effects
bool ExpressionDataSimulator::CreateDiffCoexpMatrixNoME(uint M, uint N, 
        double meanExpression, mat A, double randSdNoise, double sdNoise, 
        vector<uint> sampleIndicesInteraction) {
    // create a random data matrix
    mat D(M, N);
    // rnorm(M * N, mean = meanExpression, sd = randSdNoise)
    mt19937 engine;  // Mersenne twister random number engine
    normal_distribution<double> normDist1(meanExpression, randSdNoise);
    D.imbue( [&]() { return normDist1(engine); } );
    
    // add co-expression
    normal_distribution<double> normDist2(0, sdNoise);
    vec already_modified(M);
    already_modified.zeros();
    already_modified[0] = 1;
    for (i=0; i< M; ++i) {
      for (j=0; j < M; ++j) {
        PP->printLOG("Considering A: row" + int2str(i) + "column" + int2str(j) + "\n");
        if ((A(i, j) == 1) && (!already_modified(j))) {
          PP->printLOG("Making row" + int2str(j) + "from row" + int2str(i) + "\n");
          rowvec offset(N);
          offset.imbue( [&]() { return normDist2(engine); } )
          D(j,) = D(i,) + offset;
          already_modified(j) = 1;
        } else {
          if (already_modified(j) == 1 && !already_modified(i)) {
            // if j is already modified, we want to modify i,
            // unless i is already modified then do nothing
            rowvec offset(N);
            offset.imbue( [&]() { return normDist2(engine); } )
            D(i, ) = D(j, ) + offset;
          }
        }
      }
    }
    
    // perturb to get differential co-expression
    double n1 = N / 2;
    uniform_real_distribution<double> runif();
    uint mGenesToPerturb = sampleIndicesInteraction.size();
    for (i=0; i < mGenesToPerturb; ++i) {
      uint geneIdxInteraction = sampleIndicesInteraction(i);
      
      rowvec g0 = D(sampleIndicesInteraction, (n1 + 1):N);
      rowvec g1 = D(sampleIndicesInteraction, 1:n1);
      
      // get the group 2 gene expression and randomly order for differential coexpression
      rowvec g2 = D(geneIdxInteraction, 1:n1);
      rowvec x = shuffle(g2);
      D(geneIdxInteraction, 1:n1) = x;
    }
    
    // return a regression ready data frame
    uint dimN = D.n_cols;
    double n1 = dimN / 2;
    double n2 = dimN / 2;
    vector<string> subIds =
      c(paste("case", 1:n1, sep = ""), paste("ctrl", 1:n2, sep = ""));
    vector<unit> phenos = c(rep(1, n1), rep(0, n2));
    simulatedData = D.t() + phenos;
            
    return true;
}

// -----------------------------------------------------------------------------
// simulating X = BS + \Gamma G + U
// S = Biological group                                                                                                                   m
// G = Batch
// U = random error
// NOTE:  If you use conf=TRUE, then you must have exactly two surrogate variables in the database #
// this function only allows for confounding in the database, not confounding in the new samples #
bool ExpressionDataSimulator::SimulateData(uint n_e, uint n_db, uint n_ns,
        std::vector<std::string> sv_db, std::vector<std::string> sv_ns,
        double sd_b, double sd_gam, double sd_u, bool conf, 
        std::string distr_db, double p_b, double p_gam, double p_ov) {
  // n.e=1000 // number of variables
  // n.db=70  // sample size in database
  // n.ns=30  // sample size in newsample
  // sv.db=c("A","B") // batches
  // sv.ns=c("A","B") // batches
  // sd.b=1
  // sd.gam=1
  // sd.u=1
  // conf=FALSE
  // distr.db=NA
  // p.b=0.3
  // p.gam=0.3
  // p.ov=0.1
  uint n = n.db + n.ns;
  // Create random error
  mat U = matrix(nrow = n.e, ncol = n, rnorm(n.e * n, sd = sd.u));
  
  // Create index for database vs. new sample #
  vector<string> ind = as.factor(c(rep("db", n.db), rep("ns", n.ns)));
  
  // Create outcome, surrogate variables #
  // Use distr option to show % overlap of outcome, surrogate variables. #
  // Note that .5 means no confounding between outcome, surrogate variables. #
  
  // biological variable (fixed at 50% for each outcome)
  vector<double> S.db = c(rep(0, round(.5 * n.db)), rep(1, n.db - round(.5 * n.db)));
  vector<double> S.ns = c(rep(0, round(.5 * n.ns)), rep(1, n.ns - round(.5 * n.ns)));
  vector<double> S = c(S.db, S.ns);
  
  double len0 = sum(S.db == 0);
  double len1 = sum(S.db == 1);
  
  if (conf == FALSE) {
    // surrogate variable (no confounding in this function)
    uint n.sv.db = length(sv.db);
    double prop.db = 1 / n.sv.db;
    
    // create surrogate variables for outcome 0 in database #
    x1 = c();
    for (i in 1:n.sv.db) {
      x1 = c(x1, rep(sv.db[i], floor(prop.db * len0)));
    }
    // If the rounding has caused a problem, randomly assign to fill out vector #
    while (length(x1) != len0) {
      x1 = c(x1, sample(sv.db, 1));
    }
    
    // surrogate variables for outcome 1 will be the same #
    // this helps control for the randomly assignment - makes sure there is no #
    // added confounding #
    
    x2 = x1;
  }
  
  if (conf == TRUE) {
    x1 = c(rep("A", round(distr.db * len0)),
            rep("B", len0 - round(distr.db * len0)));
    
    x2 = c(rep("A", round((1 - distr.db) * len1)),
            rep("B", len1 - round((1 - distr.db) * len1)));
  }
  
  // create surrogate variables for outcome 0 in new samples #
  n.sv.ns = length(sv.ns);
  prop.ns = 1 / n.sv.ns;
  
  len0 = sum(S.ns == 0);
  len1 = sum(S.ns == 1);
  
  x3 = c();
  for (i in 1:n.sv.ns) {
    x3 = c(x3, rep(sv.ns[i], floor(prop.ns * len0)));
  }
  // If the rounding has caused a problem, randomly assign to fill out vector #
  while (length(x3) != len0) {
    x3 = c(x3, sample(sv.ns, 1));
  }
  
  // surrogate variables for outcome 1 will be the same #
  // this helps control for the randomly assignment - makes sure there is no #
  // added confounding #
  
  x4 = x3;
  
  G = c(x1, x2, x3, x4);
  G = t(model.matrix( ~ as.factor(G)))[-1, ];
  if (is.null(dim(G))) {
    G = matrix(G, nrow = 1, ncol = n);
  }
  
  // Determine which probes are affected by what: #
  // 30% for biological, 30% for surrogate, 10% overlap #
  // First 30% of probes will be affected by biological signal #
  ind.B = rep(0, n.e);
  ind.B[1:round(p.b * n.e)] = 1;
  // Probes 20% thru 50% will be affected by surrogate variable
  ind.Gam = rep(0, n.e);
  ind.Gam[round((p.b - p.ov) * n.e):round((p.b - p.ov + p.gam) * n.e)] = 1;
  
  // figure out dimensions for Gamma #
  
  // create parameters for signal, noise #
  B = matrix(nrow = n.e,
              ncol = 1,
              rnorm(n.e, mean = 0, sd = sd.b) * ind.B);
  Gam =
    matrix(
      nrow = n.e,
      ncol = dim(G)[1],
      rnorm(n.e * dim(G)[1], mean = 0, sd = sd.gam) * ind.Gam
    );
  
  // simulate the data #
  sim.dat = B %*% S + Gam %*% G + U;
  sim.dat = sim.dat + abs(min(sim.dat)) + 0.0001;
  
  // simulate data without batch effects #
  sim.dat.nobatch = B %*% S + U;
  sim.dat.nobatch = sim.dat.nobatch + abs(min(sim.dat)) + 0.0001;
  
  // divide parts into database, new samples #
  db = list();
  db$dat = sim.dat[, ind == "db"];
  db$datnobatch = sim.dat.nobatch[, ind == "db"];
  db$U = U[, ind == "db"];
  db$B = B;
  db$S = S[ind == "db"];
  db$Gam = Gam;
  db$G = G[ind == "db"];
  
  // Run sva on the database (will be needed for both versions of fsva) #
  //   library(sva)
  //   db$mod=model.matrix(~as.factor(db$S))
  //   db$sv=sva(db$dat,db$mod)
  new = list();
  //   new$dat = sim.dat[,ind=="ns"]
  //   new$datnobatch = sim.dat.nobatch[,ind=="ns"]
  //   new$U = U[,ind=="ns"]
  //   new$B = B
  //   new$S = S[ind=="ns"]
  //   new$Gam = Gam
  //   new$G = G[ind=="ns"]
  
  vars = list(
    n.e = n.e,
    n.db = n.db,
    n.ns = n.ns,
    sv.db = sv.db,
    sv.ns = sv.ns,
    sd.b = sd.b,
    sd.gam = sd.gam,
    sd.u = sd.u,
    conf = conf,
    distr.db = distr.db,
    p.b = p.b,
    p.gam = p.gam,
    p.ov = p.ov
  );
  
  return(list(db = db, new = new, vars = vars));
}

// -----------------------------------------------------------------------------
//###########################################################
// data simulation:
//###########################################################
bool ExpressionDataSimulator::Simulate(uint n, uint d, double pb, 
        double bias, stringtype) {
  // n = 150
  // d = 1000
  // pb = 0.05 // percentage of signal in data, 0.1 = 10% = 50 signal variables
  // type = "sva"
  // type = "pri"
  nbias = pb * d
  signal.names = paste("gene", sprintf("%04d", 1:nbias), sep = "")
  
  if (type == "sva") {
    // new simulation:
    // sd.b sort of determines how large the signals are
    // p.b=0.1 makes 10% of the variables signal
    //   bias = 0.5
    my.sim.data =
      sim_dat(
        n.e = d - 1,
        n.db = 3 * n,
        sd.b = bias,
        p.b = pb
      )$db
    data = cbind(t(my.sim.data$datnobatch), my.sim.data$S)
  } else if (type == "pri") {
    // old simulation:
    data = rnorm(n * d * 3, 0, 1)
    data = matrix(data, n * 3, d)
    data[, d] = sign(data[, d])
    signal = 1
    if (signal == 1)
      data[, 1:nbias] = data[, 1:nbias] + bias * data[, d]
  } else if (type == "inte") {
    // interaction simulation: scale-free
    g = barabasi.game(d - 1, directed = F)
    //   foo = printIGraphStats(g)
    A = get.adjacency(g)
    degrees = rowSums(A)
    myA = as.matrix(A)
    data = createDiffCoexpMatrixNoME(
      d - 1,
      n * 3,
      meanExpression = 7,
      A = myA,
      randSdNoise = 1,
      sdNoise = bias,
      1:(nbias)
    )$regressionData
  } else if (type == "er") {
    p = 0.1
    g = erdos.renyi.game(d - 1, p)
    //   foo = printIGraphStats(g)
    A = get.adjacency(g)
    // degrees = rowSums(A)
    myA = as.matrix(A)
    data = createDiffCoexpMatrixNoME(
      d - 1,
      n * 3,
      meanExpression = 7,
      A = myA,
      randSdNoise = 1,
      sdNoise = bias,
      1:(nbias)
    )$regressionData
  }
  shortname = paste(round(bias, digits = 1), pb, d, n, sep = "_")
  myfile =
    paste("data/", type, "_", shortname, "_data.Rdata", sep = "")
  
  ind = sample(3, n * 3, replace = T)
  data = data.frame(data)
  data[, d] = factor(data[, d])
  levels(data[, d]) = c(-1, 1)
  colnames(data) =
    c(paste("gene", sprintf("%04d", 1:(d - 1)), sep = ""), "pheno")
  X_train = data[ind == 1, ]
  X_holdo = data[ind == 2, ]
  X_test  = data[ind == 3, ]
  
  save(n,
       d,
       pb,
       X_train,
       X_holdo,
       X_test,
       signal.names,
       bias,
       type,
       shortname,
       file = myfile)
}
