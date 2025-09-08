import Mathlib.Tactic

open Nat

partial def smallest_factor (n : ℕ) : Option ℕ :=
  if n ≤ 1 then none
  else
    let bound := Nat.sqrt n + 1
    let rec aux (k : ℕ) : Option ℕ :=
      if k > bound then none
      else if k * k > n then none
      else if n % k = 0 then some k
      else aux (k+1)
    aux 2

partial def factor_recursive : ℕ → List ℕ
| 0 => []
| 1 => []
| n =>
  match smallest_factor n with
  | none   => [n]
  | some d =>
    let q := n / d
    d :: factor_recursive q

#eval factor_recursive 2255
