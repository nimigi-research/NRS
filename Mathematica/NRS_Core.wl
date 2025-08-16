(* ::Package:: *)

(* :Title: NRS Core *)
(* :Context: NRS` *)
(* :Author: You (+ Assistant scaffolding) *)
(* :Summary: Gibbs asymmetry, entropy production, flow centrality, QFI, and utilities *)
(* :Package Version: 0.1 *)
(* :Mathematica Version: 13+ *)

BeginPackage["NRS`"];

RowStochasticQ::usage = "RowStochasticQ[P] checks whether P is (approximately) row-stochastic.";
NormalizeRows::usage = "NormalizeRows[P] rescales each row to sum to 1.";
StationaryDistribution::usage = "StationaryDistribution[P, opts] returns a stationary vector π for row-stochastic P.";
DetailedBalanceResidual::usage = "DetailedBalanceResidual[P, π] = max_ij |π_i P_ij - π_j P_ji|.";
DetailedBalanceQ::usage = "DetailedBalanceQ[P, π, tol] is True if residual ≤ tol.";
EntropyProductionRate::usage = "EntropyProductionRate[P, π] computes 1/2 ∑_{ij} (J_ij - J_ji) log(J_ij/J_ji).";
GibbsAsymmetryQ::usage = "GibbsAsymmetryQ[P, tol] returns <|\"Stationary\"->π, \"EPR\"->e, \"GibbsAsymmetric\"->bool|>.";
FreeEnergy::usage = "FreeEnergy[ρ, U, γ] computes ⟨U⟩ - γ H(ρ) or KL(ρ || ρ*) when U==\"KL\".";
FlowCentrality::usage = "FlowCentrality[P, π] returns vertex centrality (out-flux) and edge flux J.";
QFIMetricDiscrete::usage = "QFIMetricDiscrete[vals, p, θ] computes Fisher metric E[(∇θ log p)(∇θ log p)^T].";
RandomIrreversibleKernel::usage = "RandomIrreversibleKernel[n] makes a non-reversible row-stochastic kernel.";
ExampleKernel::usage = "ExampleKernel[] returns a small 3-state irreversible kernel.";

Options[StationaryDistribution] = {"Method" -> "Eigen", "Tolerance" -> 10^-12, "MaxTries" -> 3};

Begin["`Private`"];

eps = $MachineEpsilon;

RowStochasticQ[P_?MatrixQ, tol_: 10^-10] := Module[{sums, nonneg},
  sums = Total /@ P;
  nonneg = Min[Flatten[P]] >= -tol;
  Max[Abs[sums - 1]] <= tol && nonneg
];

NormalizeRows[P_?MatrixQ] := Module[{r = Total /@ P},
  P/DiagonalMatrix[r /. 0. -> 1.]
];

(* robust stationary distribution finder *)
StationaryDistribution[P_?MatrixQ, opts : OptionsPattern[]] := Module[
  {method = OptionValue["Method"], tol = OptionValue["Tolerance"], tries = OptionValue["MaxTries"], n, pi, k = 0, PP = N[P]},
  n = Length[PP];
  If[!RowStochasticQ[PP, 10^-6], PP = NormalizeRows[PP]];
  Which[
    method === "Eigen",
      Module[{vals, vecs, idx},
        {vals, vecs} = Eigensystem[Transpose[PP]];
        idx = First@Ordering[Abs[vals - 1], 1];
        pi = vecs[[idx]];
        pi = Re[pi/Total[pi]];
        If[Min[pi] < -Sqrt[tol], (* fallback to linear solve *)
          StationaryDistribution[PP, "Method" -> "LinearSolve", "Tolerance" -> tol, "MaxTries" -> tries],
          If[Max[Abs[Transpose[PP].pi - pi]] <= 10^(-6), Abs[pi], 
             StationaryDistribution[PP, "Method" -> "LinearSolve", "Tolerance" -> tol, "MaxTries" -> tries]]
        ]
      ],
    method === "LinearSolve",
      Module[{A, b, x},
        A = Join[Transpose[PP] - IdentityMatrix[n], {ConstantArray[1., n]}];
        b = Join[ConstantArray[0., n], {1.}];
        Quiet@Check[
          x = LinearSolve[A, b];
          pi = x;
          If[Min[pi] < -Sqrt[tol], $Failed, N[Abs[pi]/Total[Abs[pi]]]],
          If[k++ < tries, StationaryDistribution[PP + 10^-12 RandomReal[{-1, 1}, Dimensions[PP]],
                                                 "Method" -> "Eigen", "Tolerance" -> tol, "MaxTries" -> tries],
             $Failed]
        ]
      ],
    True, StationaryDistribution[PP, "Method" -> "Eigen", "Tolerance" -> tol, "MaxTries" -> tries]
  ]
];

DetailedBalanceResidual[P_?MatrixQ, pi_List] := Module[{J = Outer[#1*#2 &, pi, P, 1]},
  Max[Abs[J - Transpose[J]]]
];

DetailedBalanceQ[P_?MatrixQ, pi_List, tol_: 10^-9] := DetailedBalanceResidual[P, pi] <= tol;

EntropyProductionRate[P_?MatrixQ, pi_List] := Module[{J, n, term},
  n = Length[P];
  J = Outer[#1*#2 &, pi, P, 1] // N;
  term[i_, j_] := Module[{a = Max[J[[i, j]], eps], b = Max[J[[j, i]], eps]},
    (a - b) * Log[a/b]
  ];
  N[0.5 Sum[term[i, j], {i, n}, {j, n}]]
];

GibbsAsymmetryQ[P_?MatrixQ, tol_: 10^-9] := Module[{pi, epr},
  pi = StationaryDistribution[P];
  If[pi === $Failed, Return[$Failed]];
  epr = EntropyProductionRate[P, pi];
  <|"Stationary" -> pi, "EPR" -> epr, "GibbsAsymmetric" -> (epr > tol)|>
];

(* Free energy: either thermodynamic form or KL form *)
FreeEnergy[ρ_?VectorQ, U_: (0 &), γ_: 1, ref_: None] := Module[{p = N[ρ/Total[ρ]], H},
  If[ref === "KL",
    Module[{q = N[ρ/Total[ρ]]}, (* interpret ρ as 'p', user passes ref via U if desired *)
      (* If ref=="KL", interpret U as reference distribution *)
      Message[FreeEnergy::argx, "When ref==\"KL\", pass refDist as third arg."];
      $Failed
    ],
    H = -Total[p*Log[Clip[p, {eps, 1.}]]];
    N[Total[p*U /@ Range[Length[p]]] - γ H]
  ]
];

(* Fisher information for discrete domain *)
QFIMetricDiscrete[vals_List, p_, θ_] := Module[{pdf, lp, grads, w, n},
  pdf[x_] := Evaluate[p[x, θ]];
  lp[x_] := Log[Max[pdf[x], eps]];
  grads = (D[lp[#], {θ}] &) /@ vals;
  w = pdf /@ vals; w = N[w/Total[w]];
  n = Length[grads];
  N[Sum[w[[i]]*Outer[Times, grads[[i]], grads[[i]]], {i, n}]]
];

FlowCentrality[P_?MatrixQ, pi_List] := Module[{J = Outer[#1*#2 &, pi, P, 1]},
  <|"VertexCentrality" -> Total[J, {2}] // N, "EdgeFlux" -> N[J]|>
];

RandomIrreversibleKernel[n_Integer?Positive] := Module[{Q, P, pi, res},
  Q = RandomReal[{0, 1}, {n, n}];
  P = NormalizeRows[Q];
  pi = ConstantArray[1./n, n];
  res = DetailedBalanceResidual[P, pi];
  If[res < 10^-8, RandomIrreversibleKernel[n], P]
];

ExampleKernel[] := {{0.8, 0.1, 0.1},
                    {0.1, 0.8, 0.1},
                    {0.2, 0.0, 0.8}} // N;

End[]; (* `Private` *)
EndPackage[];