/-
  NRS — Kernel & EPR core (Lean 4 + mathlib)
  Finite-state stochastic kernels, stationary distributions, detailed balance, EPR.
-/

import Mathlib.Data.Matrix.Basic
import Mathlib.Analysis.SpecialFunctions.Log
import Mathlib.Tactic
import Mathlib.Data.Finset.Basic
import Mathlib.Algebra.BigOperators.Ring

open BigOperators Matrix

namespace NRS

variable {X : Type*} [Fintype X] [DecidableEq X]

abbrev ℝ≥ := Real

/-- A row-stochastic kernel on a finite type. -/
structure StochMatrix (X : Type*) [Fintype X] [DecidableEq X] where
  P      : Matrix X X ℝ≥
  nonneg : ∀ i j, 0 ≤ P i j
  rowsum : ∀ i, (∑ j, P i j) = 1

attribute [simp] StochMatrix.rowsum

/-- A stationary distribution for `P`. -/
structure Stationary (P : StochMatrix X) where
  π       : X → ℝ≥
  nonneg  : ∀ i, 0 ≤ π i
  sum_one : (∑ i, π i) = 1
  stat    : (fun j => ∑ i, (π i) * P.P i j) = π

/-- Detailed balance w.r.t. π. -/
def detailedBalance {P : StochMatrix X} (π : X → ℝ≥) : Prop :=
  ∀ i j, (π i) * P.P i j = (π j) * P.P j i

/-- Edge flux J_{ij} = π_i P_{ij}. -/
def flux {P : StochMatrix X} (π : X → ℝ≥) (i j : X) : ℝ≥ := (π i) * P.P i j

/-- Entropy production rate (finite-state). -/
noncomputable def EPR {P : StochMatrix X} (S : Stationary P) : ℝ≥ :=
  (1/2 : ℝ≥) *
    ∑ i, ∑ j,
      let a := flux (S.π) i j
      let b := flux (S.π) j i
      (a - b) * Real.log (a / b)

-- Positivity assumptions (strictly positive flux avoids log 0).
def positiveFlux {P : StochMatrix X} (S : Stationary P) : Prop :=
  ∀ i j, 0 < flux S.π i j

/-- EPR ≥ 0 (skeleton; proof via log-sum inequality / KL symmetry). -/
theorem epr_nonneg {P : StochMatrix X} (S : Stationary P)
    (pos : positiveFlux S) : 0 ≤ EPR (P:=P) S := by
  -- Idea: EPR = 1/2 ( KL(J || Jᵗ) + KL(Jᵗ || J) ) where J_{ij} = π_i P_{ij}.
  -- Each KL ≥ 0.
  -- A full proof rewrites the double sum and applies convexity of x ↦ x log x.
  sorry

/-- EPR = 0 iff detailed balance (skeleton). -/
theorem epr_eq_zero_iff_db {P : StochMatrix X} (S : Stationary P)
    (pos : positiveFlux S) :
    EPR (P:=P) S = 0 ↔ detailedBalance (P:=P) S.π := by
  -- `→`: zero symmetric KL ⇒ J = Jᵗ, i.e. π_i P_{ij} = π_j P_{ji}.
  -- `←`: DB makes each summand vanish.
  sorry

end NRS