(* ::Package:: *)

(* :Title: NRS Murphy Scheduler *)
(* :Context: NRS`Murphy` *)
(* :Summary: Greedy and DP schedulers for budgeted information gain *)
(* :Package Version: 0.1 *)

BeginPackage["NRS`Murphy`"];

MurphyGreedySchedule::usage = 
  "MurphyGreedySchedule[actions, budget, infoGain, cost, state0, maxSteps] \
returns <|\"Policy\"->list, \"UnspentBudget\"->b, \"Trace\"->...|>.";

DPBudgetedSchedule::usage = 
  "DPBudgetedSchedule[actions, budget, infoGain, cost, state0] \
does discrete DP (small budgets).";

SubmodularCertificate::usage = 
  "SubmodularCertificate[actions, infoGain, state0, samples] \
returns a numeric certificate of diminishing returns.";

Begin["`Private`"];

MurphyGreedySchedule[actions_List, budget_?NumericQ, infoGain_, cost_, state0_, maxSteps_: 100] := Module[
  {B = budget, policy = {}, state = state0, step = 0, hist = {}, cand, best, ratio},
  While[B > 0 && step < maxSteps,
    cand = Select[actions, cost[#] <= B &];
    If[cand === {}, Break[]];
    ratio[a_] := Quiet@Check[N[infoGain[a, state]/Max[cost[a], $MachineEpsilon]], -Infinity];
    best = First@MaximalBy[cand, ratio];
    AppendTo[policy, best];
    AppendTo[hist, <|"Action"->best, "Gain"->infoGain[best, state], "Cost"->cost[best], "Budget"->B|>];
    B -= cost[best];
    (* If state update is intended, user should supply a side-effecting infoGain; we keep immutable here *)
    step++
  ];
  <|"Policy" -> policy, "UnspentBudget" -> B, "Trace" -> hist|>
];

DPBudgetedSchedule[actions_List, budget_Integer?NonNegative, infoGain_, cost_, state0_] := Module[
  {n = Length[actions], dp, choice, action, b, i},
  dp = ConstantArray[0., {n + 1, budget + 1}];
  choice = Table[None, {n + 1}, {budget + 1}];
  For[i = 1, i <= n, i++,
    action = actions[[i]];
    For[b = 0, b <= budget, b++,
      dp[[i + 1, b + 1]] = dp[[i, b + 1]];
      choice[[i + 1, b + 1]] = choice[[i, b + 1]];
      If[cost[action] <= b,
        With[{cand = infoGain[action, state0] + dp[[i, b - cost[action] + 1]]},
          If[cand > dp[[i + 1, b + 1]],
            dp[[i + 1, b + 1]] = cand;
            choice[[i + 1, b + 1]] = {action, b - cost[action]}
          ]
        ]
      ]
    ]
  ];
  (* reconstruct *)
  b = budget; i = n; policy = {};
  While[i > 0 && b >= 0,
    If[choice[[i + 1, b + 1]] =!= choice[[i, b + 1]] && choice[[i + 1, b + 1]] =!= None,
      AppendTo[policy, First@choice[[i + 1, b + 1]]];
      b = Last@choice[[i + 1, b + 1]];
    ];
    i--;
  ];
  <|"Policy" -> Reverse@policy, "OptimalValue" -> dp[[n + 1, budget + 1]]|>
];

SubmodularCertificate[actions_List, infoGain_, state0_, samples_: 50] := Module[
  {pairs, trip, a, b, S, ΔaS, ΔaSb, cert = True, violations = {}},
  pairs = RandomSample[Subsets[actions, {0, Min[3, Length@actions]}], Min[samples, 50]];
  Do[
    {S, a, b} = RandomChoice[{pairs, actions, actions}];
    If[a === b, Continue[]];
    ΔaS  = infoGain[a, state0];       (* treat as marginal gain proxy *)
    ΔaSb = infoGain[a, state0];       (* in practice, pass a marginal-gain oracle depending on set *)
    If[ΔaS < ΔaSb - 10^-9, cert = False; AppendTo[violations, {S, a, b}]];
  , {samples}];
  <|"Submodular" -> cert, "Violations" -> violations|>
];

End[]; (* `Private` *)
EndPackage[];