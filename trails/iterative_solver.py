from __future__ import annotations

from typing import Literal, Optional
import warnings

import numpy as np
from scipy import sparse as sp
from scipy.sparse.linalg import LinearOperator, gmres, spilu


def _build_jacobi_preconditioner(
    A_csc: sp.csc_matrix,
) -> Optional[LinearOperator]:
    """Build a Jacobi preconditioner ``D^-1`` for sparse matrix ``A``."""
    diagonal = np.asarray(A_csc.diagonal(), dtype=np.float64)
    if np.any(diagonal == 0.0):
        return None

    inv_diag = 1.0 / diagonal
    return LinearOperator(
        shape=A_csc.shape,
        matvec=lambda x: inv_diag * x,
        dtype=np.float64,
    )


def _build_ilu_preconditioner(
    A_csc: sp.csc_matrix,
    *,
    drop_tol: float,
    fill_factor: float,
) -> Optional[LinearOperator]:
    """Build an ILU preconditioner, returning ``None`` if factorization fails."""
    try:
        ilu = spilu(A_csc, drop_tol=float(drop_tol), fill_factor=float(fill_factor))
    except Exception:
        return None

    return LinearOperator(
        shape=A_csc.shape,
        matvec=ilu.solve,
        dtype=np.float64,
    )


def solve_many_rhs_jacobi_gmres(
    A_csc: sp.csc_matrix,
    B: np.ndarray,
    *,
    rtol: float = 1e-6,
    atol: float = 0.0,
    restart: int | None = 50,
    maxiter: int | None = 300,
    use_guess: bool = True,
    preconditioner_mode: Literal["jacobi", "ilu", "none"] = "jacobi",
    ilu_drop_tol: float = 1e-4,
    ilu_fill_factor: float = 10.0,
) -> np.ndarray:
    """Solve ``A X = B`` with GMRES and configurable preconditioning."""
    if not sp.isspmatrix_csc(A_csc):
        A_csc = A_csc.tocsc()

    if A_csc.dtype != np.float64:
        A_csc = A_csc.astype(np.float64)

    # Keep matrix in canonical CSC form for stable/fast iterative solves.
    A_csc.sum_duplicates()
    A_csc.eliminate_zeros()
    A_csc.sort_indices()

    B = np.asarray(B)
    if B.ndim != 2:
        raise ValueError("B must be 2D (n, k)")
    if B.dtype != np.float64:
        B = B.astype(np.float64, copy=False)

    n, k = B.shape
    if A_csc.shape != (n, n):
        raise ValueError(f"Shape mismatch: A {A_csc.shape}, B {B.shape}")

    if preconditioner_mode not in {"jacobi", "ilu", "none"}:
        raise ValueError(
            "preconditioner_mode must be one of {'jacobi', 'ilu', 'none'}"
        )

    preconditioner: Optional[LinearOperator]
    if preconditioner_mode == "none":
        preconditioner = None
    elif preconditioner_mode == "jacobi":
        preconditioner = _build_jacobi_preconditioner(A_csc)
    else:
        preconditioner = _build_ilu_preconditioner(
            A_csc,
            drop_tol=float(ilu_drop_tol),
            fill_factor=float(ilu_fill_factor),
        )
        if preconditioner is None:
            warnings.warn(
                "ILU preconditioner failed; falling back to Jacobi preconditioner.",
                RuntimeWarning,
                stacklevel=2,
            )
            preconditioner = _build_jacobi_preconditioner(A_csc)

    X = np.empty((n, k), dtype=np.float64)
    guess: np.ndarray | None = None

    for j in range(k):
        rhs = B[:, j]
        x0 = guess if (use_guess and guess is not None) else None

        try:
            sol, info = gmres(
                A_csc,
                rhs,
                x0=x0,
                rtol=rtol,
                atol=atol,
                restart=restart,
                maxiter=maxiter,
                M=preconditioner,
            )
        except TypeError:
            # Compatibility fallback for SciPy versions with `tol` instead of `rtol`.
            sol, info = gmres(
                A_csc,
                rhs,
                x0=x0,
                tol=rtol,
                atol=atol,
                restart=restart,
                maxiter=maxiter,
                M=preconditioner,
            )

        if info != 0:
            raise RuntimeError(
                "GMRES failed to converge "
                f"(rhs_col={j}, info={info}, rtol={rtol}, maxiter={maxiter})"
            )

        sol_arr = np.asarray(sol, dtype=np.float64)
        if not sol_arr.shape:
            sol_arr = sol_arr.reshape((1,))
        X[:, j] = sol_arr

        if use_guess:
            guess = sol_arr

    return X
