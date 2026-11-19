"""
Axis code generator for spatial zipcodes.

Implements a randomized greedy algorithm with backtracking to generate
1D sequences of binary codes that satisfy distance constraints.
"""

import logging
import random
from typing import List, Dict, Tuple, Optional
import numpy as np

from .sequence_qc import (
    calculate_g_content,
    find_max_homopolymer,
    max_right_edge_by_base,
    terminal_run,
)

logger = logging.getLogger(__name__)


def hamming_distance(code1: List[int], code2: List[int]) -> int:
    """
    Calculate Hamming distance between two binary codes.
    
    Args:
        code1: First binary code (list of 0s and 1s)
        code2: Second binary code (list of 0s and 1s)
        
    Returns:
        Hamming distance (number of differing bits)
    """
    if len(code1) != len(code2):
        raise ValueError(f"Codes must have same length: {len(code1)} vs {len(code2)}")
    return sum(b1 != b2 for b1, b2 in zip(code1, code2))


def random_binary_code(length: int, rng: random.Random) -> List[int]:
    """
    Generate a random binary code of specified length.
    
    Args:
        length: Number of bits
        rng: Random number generator
        
    Returns:
        List of 0s and 1s
    """
    return [rng.randint(0, 1) for _ in range(length)]


def flip_neighbour_code(previous: List[int], n_flips: int, rng: random.Random) -> List[int]:
    """
    Copy a binary code and flip a number of distinct bits.

    Args:
        previous: Code to start from (left unchanged)
        n_flips: Number of distinct bits to flip
        rng: Random number generator

    Returns:
        New list of 0s and 1s at Hamming distance n_flips from previous
    """
    if not 0 <= n_flips <= len(previous):
        raise ValueError(f"n_flips must be in [0, {len(previous)}], got {n_flips}")
    candidate = list(previous)
    for bit in rng.sample(range(len(previous)), n_flips):
        candidate[bit] ^= 1
    return candidate


CANDIDATE_MODES = ("random", "flip")


class AxisCodeGenerator:
    """
    Generates 1D axis codes with specified distance constraints.
    
    Uses a randomized greedy algorithm with backtracking to construct
    codes that satisfy hard constraints (far-features minimum distance)
    and approximate soft constraints (target distances for nearby positions).

    Candidate modes:
        "random": every candidate is a fresh random code. At 36 bits a random
            code is about 18 bits from any other, so small targets are missed.
        "flip": every candidate is the preceding code with a few bits flipped
            (counts drawn from flip_counts). The first candidate that meets all
            targets exactly is taken; otherwise the best-scoring one.
    """
    
    def __init__(
        self,
        n_positions: int,
        bit_length: int,
        target_distances: List[int],
        far_features_min_distance: int,
        max_attempts_per_position: int = 10000,
        backtrack_depth: int = 5,
        max_backtracks: int = 100,
        tolerance: float = 0.5,
        random_seed: Optional[int] = None,
        encoder=None,
        g_content_max: Optional[float] = None,
        max_homopolymer_length: Optional[int] = None,
        spacer: str = "",
        junction_role: Optional[str] = None,
        partner_right_edge: Optional[Dict[str, int]] = None,
        downstream_start_runs: Optional[Dict[str, int]] = None,
        candidate_mode: str = "random",
        flip_counts: Tuple[int, ...] = (1, 2, 3),
        require_tolerance: bool = False,
    ):
        """
        Initialize the axis code generator.
        
        Args:
            n_positions: Number of positions along the axis
            bit_length: Length of binary codes
            target_distances: Target Hamming distances for position differences 1,2,3,4
            far_features_min_distance: Minimum distance for positions >= 5 apart
            max_attempts_per_position: Max random samples per position before backtracking
            backtrack_depth: Number of positions to backtrack when stuck
            max_backtracks: Maximum number of backtrack events before giving up
            tolerance: Tolerance for soft constraint violations (fraction)
            random_seed: Random seed for reproducibility
            encoder: Bit-to-DNA encoder. Required when a G cap or junction
                budget is set. Composition is judged on the encoded sequence.
            g_content_max: Maximum G fraction of the encoded axis sequence.
                Candidates above this are rejected before any distance is scored.
            max_homopolymer_length: Maximum junction homopolymer, in bases,
                allowed on the concatenated barcode.
            spacer: Bases inserted between the column code and the row code.
            junction_role: "left" for the column code, "right" for the row code.
                Left codes keep a free starting base for the other axis.
                Right codes are rejected when their leading run, added to the
                matching column-code edge, would exceed max_homopolymer_length.
            partner_right_edge: For role "right", longest right-edge run of each
                base across the already chosen column codes (spacer included).
            downstream_start_runs: Longest leading run of each base among sequences
                that will be appended after this code (e.g. HSDB bar_2). The
                code's trailing run plus that start must fit in
                max_homopolymer_length.
            candidate_mode: "random" (fresh random codes) or "flip" (bit flips
                of the preceding code). Position 0 is always random.
            flip_counts: Numbers of bits a "flip" candidate may flip, drawn
                uniformly per candidate.
            require_tolerance: If True, a position with no candidate inside the
                tolerance triggers backtracking. If False, the best candidate
                found is kept even when it misses the tolerance.
        """
        if candidate_mode not in CANDIDATE_MODES:
            raise ValueError(f"candidate_mode must be one of {CANDIDATE_MODES}, got {candidate_mode!r}")
        if candidate_mode == "flip":
            if not flip_counts or any(not 1 <= k <= bit_length for k in flip_counts):
                raise ValueError(f"flip_counts must be non-empty values in [1, {bit_length}], got {flip_counts}")
        if g_content_max is not None and encoder is None:
            raise ValueError("encoder is required when g_content_max is set")
        if junction_role is not None and junction_role not in ("left", "right"):
            raise ValueError(f"junction_role must be 'left' or 'right', got {junction_role!r}")
        if junction_role is not None:
            if encoder is None:
                raise ValueError("encoder is required when junction_role is set")
            if max_homopolymer_length is None or max_homopolymer_length < 1:
                raise ValueError(
                    "max_homopolymer_length must be a positive integer when "
                    "junction_role is set"
                )
        if junction_role == "right" and partner_right_edge is None:
            raise ValueError("partner_right_edge is required for junction_role 'right'")
        if junction_role == "left" and not (spacer or "") and max_homopolymer_length < 2:
            raise ValueError(
                "max_homopolymer_length must be at least 2 when the spacer is empty, "
                "so the column code can leave room for the row code"
            )
        if downstream_start_runs is not None:
            if encoder is None:
                raise ValueError("encoder is required when downstream_start_runs is set")
            if max_homopolymer_length is None or max_homopolymer_length < 1:
                raise ValueError(
                    "max_homopolymer_length must be a positive integer when "
                    "downstream_start_runs is set"
                )
        if g_content_max is not None and not 0 <= g_content_max <= 1:
            raise ValueError(f"g_content_max must be in [0, 1], got {g_content_max}")

        self.n_positions = n_positions
        self.bit_length = bit_length
        self.target_distances = target_distances
        self.far_features_min_distance = far_features_min_distance
        self.max_attempts = max_attempts_per_position
        self.backtrack_depth = backtrack_depth
        self.max_backtracks = max_backtracks
        self.tolerance = tolerance
        self.encoder = encoder
        self.g_content_max = g_content_max
        self.max_homopolymer_length = max_homopolymer_length
        self.spacer = spacer or ""
        self.junction_role = junction_role
        self.partner_right_edge = partner_right_edge or {}
        self.downstream_start_runs = dict(downstream_start_runs or {})
        self.candidate_mode = candidate_mode
        self.flip_counts = tuple(flip_counts)
        self.require_tolerance = require_tolerance
        
        self.rng = random.Random(random_seed)
        self.codes: List[List[int]] = []
        self._right_edges: List[Tuple[str, int]] = []
        self._blocked_counts = {base: 0 for base in "ACGT"}
        self.rejected_for_g = 0
        self.rejected_for_homopolymer = 0
        self.rejected_for_junction = 0
        self.rejected_for_downstream = 0
        
        logger.info(
            f"Initialized AxisCodeGenerator: n={n_positions}, L={bit_length}, "
            f"targets={target_distances}, far_min={far_features_min_distance}, "
            f"mode={candidate_mode}"
        )
    
    def _target_distance_for_diff(self, diff: int) -> int:
        """Get target distance for a position difference."""
        if diff == 0:
            return 0
        elif 1 <= diff <= 4:
            return self.target_distances[diff - 1]
        else:  # diff >= 5
            return self.far_features_min_distance
    
    def _check_hard_constraints(self, candidate: List[int], position: int) -> bool:
        """
        Check if candidate satisfies hard constraints.
        
        Hard constraint: minimum distance to all far-away positions (diff >= 5).
        
        Args:
            candidate: Candidate binary code
            position: Current position index
            
        Returns:
            True if all hard constraints satisfied
        """
        for prev_pos in range(position):
            diff = position - prev_pos
            if diff >= 5:
                dist = hamming_distance(candidate, self.codes[prev_pos])
                if dist < self.far_features_min_distance:
                    return False
        return True
    
    def _score_soft_constraints(self, candidate: List[int], position: int) -> float:
        """
        Score candidate based on soft constraints.
        
        Soft constraints: target distances for nearby positions (diff 1-4).
        Lower score is better.
        
        Args:
            candidate: Candidate binary code
            position: Current position index
            
        Returns:
            Score (sum of squared errors from target distances)
        """
        score = 0.0
        for prev_pos in range(position):
            diff = position - prev_pos
            if 1 <= diff <= 4:
                target = self.target_distances[diff - 1]
                actual = hamming_distance(candidate, self.codes[prev_pos])
                error = actual - target
                score += error ** 2
        return score
    
    def _is_acceptable(self, candidate: List[int], position: int) -> bool:
        """
        Check if candidate is acceptable (satisfies hard + soft constraints).
        
        Args:
            candidate: Candidate binary code
            position: Current position index
            
        Returns:
            True if acceptable
        """
        # Must satisfy hard constraints
        if not self._check_hard_constraints(candidate, position):
            return False
        
        # Check soft constraints within tolerance
        score = self._score_soft_constraints(candidate, position)
        
        # Calculate acceptable threshold based on position and tolerance
        # More lenient for positions with more constraints
        n_soft_constraints = min(position, 4)
        if n_soft_constraints == 0:
            return True
        
        # Average squared error per constraint
        avg_sq_error = score / n_soft_constraints
        
        # Accept if average error is within tolerance
        max_acceptable_error = (self.tolerance * max(self.target_distances)) ** 2
        
        return avg_sq_error <= max_acceptable_error

    def _passes_composition(self, candidate: List[int]) -> bool:
        """
        Reject a candidate whose DNA breaks a composition rule.

        The G check runs before any Hamming distance is computed. Junction
        budgets run next, still before distance scoring. The junction run
        itself is measured later on the concatenated barcode; these budgets
        keep that measurement inside the limit.
        """
        if self.encoder is None:
            return True

        dna = self.encoder.encode(candidate)
        if self.g_content_max is not None:
            if calculate_g_content(dna) > self.g_content_max:
                self.rejected_for_g += 1
                return False

        if self.max_homopolymer_length is not None:
            homo_len, _ = find_max_homopolymer(dna)
            if homo_len > self.max_homopolymer_length:
                self.rejected_for_homopolymer += 1
                return False

        if self.junction_role == "left":
            base, length = terminal_run(dna + self.spacer, at_end=True)
            limit = self.max_homopolymer_length
            # An empty spacer has to leave one base of room, so a row code can
            # start with that base and the junction run still fits. A spacer
            # that already spends the whole budget cannot leave that room.
            if length > self._max_left_edge() or (
                length >= limit and self._blocks_last_start_base(base)
            ):
                self.rejected_for_junction += 1
                return False
        elif self.junction_role == "right":
            base, length = terminal_run(dna, at_end=False)
            edge = self.partner_right_edge.get(base, 0)
            if edge > 0 and edge + length > self.max_homopolymer_length:
                self.rejected_for_junction += 1
                return False

        if self.downstream_start_runs:
            base, length = terminal_run(dna, at_end=True)
            start = self.downstream_start_runs.get(base, 0)
            if start > 0 and length + start > self.max_homopolymer_length:
                self.rejected_for_downstream += 1
                return False
        return True

    def _max_left_edge(self) -> int:
        """Longest right-edge run a column code may contribute to the junction."""
        limit = self.max_homopolymer_length
        if not self.spacer:
            return limit - 1
        _, spacer_run = terminal_run(self.spacer, at_end=True)
        if spacer_run >= limit:
            return limit
        return limit - 1

    def _blocks_last_start_base(self, base: str) -> bool:
        """True when accepting this edge would forbid every row-code start base."""
        blocked = {b for b, count in self._blocked_counts.items() if count > 0}
        return base not in blocked and len(blocked) == 3

    def _remember_left_edge(self, candidate: List[int]) -> None:
        if self.junction_role != "left":
            return
        base, length = terminal_run(self.encoder.encode(candidate) + self.spacer, at_end=True)
        self._right_edges.append((base, length))
        if length >= self.max_homopolymer_length:
            self._blocked_counts[base] += 1

    def _forget_left_edges(self, count: int) -> None:
        if self.junction_role != "left" or count <= 0:
            return
        for base, length in self._right_edges[-count:]:
            if length >= self.max_homopolymer_length:
                self._blocked_counts[base] -= 1
        del self._right_edges[-count:]

    def _generate_position(self, position: int) -> Optional[List[int]]:
        """
        Generate code for a single position.
        
        Args:
            position: Position index
            
        Returns:
            Binary code if successful, None if failed
        """
        best_candidate = None
        best_score = float('inf')
        flip = self.candidate_mode == "flip" and position > 0
        
        for attempt in range(self.max_attempts):
            candidate = self._draw_candidate(position)

            # Composition before distances: a high-G code is never scored.
            if not self._passes_composition(candidate):
                continue

            if self._check_hard_constraints(candidate, position):
                score = self._score_soft_constraints(candidate, position)
                
                if score < best_score:
                    best_candidate = candidate
                    best_score = score
                
                # Flip candidates hit the targets exactly often enough to wait
                # for a perfect one; the loose tolerance would accept drift.
                if flip:
                    if score == 0:
                        return candidate
                elif self._is_acceptable(candidate, position):
                    return candidate
        
        if best_candidate is None:
            return None
        if not self._is_acceptable(best_candidate, position):
            if self.require_tolerance:
                logger.debug(
                    f"Position {position}: no candidate within tolerance "
                    f"(best score {best_score:.2f}); backtracking"
                )
                return None
            logger.debug(
                f"Position {position}: accepting best candidate with score {best_score:.2f} "
                f"after {self.max_attempts} attempts"
            )
        return best_candidate

    def _draw_candidate(self, position: int) -> List[int]:
        """Draw one candidate code for a position according to candidate_mode."""
        if self.candidate_mode == "flip" and position > 0:
            n_flips = self.rng.choice(self.flip_counts)
            return flip_neighbour_code(self.codes[position - 1], n_flips, self.rng)
        return random_binary_code(self.bit_length, self.rng)
    
    def generate(self) -> List[List[int]]:
        """
        Generate axis codes for all positions.
        
        Uses greedy generation with backtracking when stuck.
        
        Returns:
            List of binary codes (one per position)
            
        Raises:
            RuntimeError: If unable to generate valid codes
        """
        self.codes = []
        self._right_edges = []
        self._blocked_counts = {base: 0 for base in "ACGT"}
        position = 0
        backtrack_count = 0
        
        logger.info(f"Starting axis code generation for {self.n_positions} positions...")
        
        while position < self.n_positions:
            if position % 10 == 0 and position > 0:
                logger.info(f"Generated codes for {position}/{self.n_positions} positions")
            
            candidate = self._generate_position(position)
            
            if candidate is not None:
                self.codes.append(candidate)
                self._remember_left_edge(candidate)
                position += 1
            else:
                # Backtrack
                backtrack_count += 1
                if backtrack_count > self.max_backtracks:
                    raise RuntimeError(
                        f"Failed to generate codes: exceeded {self.max_backtracks} backtracks"
                    )
                
                backtrack_amount = min(self.backtrack_depth, position)
                if backtrack_amount == 0:
                    raise RuntimeError(
                        f"Failed to generate code for position 0 "
                        f"(bit_length={self.bit_length} may be too small)"
                    )
                
                logger.warning(
                    f"Backtracking {backtrack_amount} positions from position {position}"
                )
                self._forget_left_edges(backtrack_amount)
                self.codes = self.codes[:-backtrack_amount]
                position -= backtrack_amount
        
        logger.info(
            f"Successfully generated {len(self.codes)} codes "
            f"(backtracks: {backtrack_count})"
        )
        
        return self.codes
    
    def verify_constraints(self) -> Dict[str, any]:
        """
        Verify generated codes satisfy constraints and compute statistics.
        
        Returns:
            Dictionary with verification results and statistics
        """
        if not self.codes:
            raise ValueError("No codes generated yet")
        
        results = {
            "n_positions": len(self.codes),
            "bit_length": self.bit_length,
            "hard_constraint_violations": 0,
            "positions_outside_tolerance": sum(
                not self._is_acceptable(self.codes[i], i) for i in range(len(self.codes))
            ),
            "soft_constraint_stats": {},
            "distance_distribution": {},
        }
        
        # Check all pairwise distances
        for i in range(len(self.codes)):
            for j in range(i + 1, len(self.codes)):
                diff = j - i
                dist = hamming_distance(self.codes[i], self.codes[j])
                
                # Record in distribution
                if diff not in results["distance_distribution"]:
                    results["distance_distribution"][diff] = []
                results["distance_distribution"][diff].append(dist)
                
                # Check hard constraint
                if diff >= 5 and dist < self.far_features_min_distance:
                    results["hard_constraint_violations"] += 1
                
                # Check soft constraints
                if 1 <= diff <= 4:
                    target = self.target_distances[diff - 1]
                    error = dist - target
                    
                    if diff not in results["soft_constraint_stats"]:
                        results["soft_constraint_stats"][diff] = {
                            "target": target,
                            "errors": [],
                        }
                    results["soft_constraint_stats"][diff]["errors"].append(error)
        
        # Compute summary statistics
        for diff, stats in results["soft_constraint_stats"].items():
            errors = stats["errors"]
            stats["mean_error"] = np.mean(errors)
            stats["std_error"] = np.std(errors)
            stats["mean_abs_error"] = np.mean(np.abs(errors))
            stats["actual_mean_distance"] = stats["target"] + stats["mean_error"]
        
        # Summary of distance distribution
        for diff, distances in results["distance_distribution"].items():
            results["distance_distribution"][diff] = {
                "min": min(distances),
                "max": max(distances),
                "mean": np.mean(distances),
                "std": np.std(distances),
            }
        
        return results


def generate_axis_codes(
    n_positions: int,
    bit_length: int,
    target_distances: List[int],
    far_features_min_distance: int,
    max_attempts_per_position: int = 10000,
    backtrack_depth: int = 5,
    max_backtracks: int = 100,
    tolerance: float = 0.5,
    random_seed: Optional[int] = None,
    encoder=None,
    g_content_max: Optional[float] = None,
    max_homopolymer_length: Optional[int] = None,
    spacer: str = "",
    junction_role: Optional[str] = None,
    partner_right_edge: Optional[Dict[str, int]] = None,
    downstream_start_runs: Optional[Dict[str, int]] = None,
    candidate_mode: str = "random",
    flip_counts: Tuple[int, ...] = (1, 2, 3),
    require_tolerance: bool = False,
) -> Tuple[List[List[int]], Dict[str, any]]:
    """
    Generate and verify axis codes.
    
    Convenience function that creates a generator, generates codes,
    and verifies constraints.
    
    Args:
        n_positions: Number of positions along the axis
        bit_length: Length of binary codes
        target_distances: Target Hamming distances for position differences 1,2,3,4
        far_features_min_distance: Minimum distance for positions >= 5 apart
        max_attempts_per_position: Max random samples per position
        backtrack_depth: Number of positions to backtrack when stuck
        max_backtracks: Maximum number of backtrack events before giving up
        tolerance: Tolerance for soft constraint violations
        random_seed: Random seed for reproducibility
        encoder: Bit-to-DNA encoder used for composition checks
        g_content_max: Maximum G fraction of each encoded axis sequence
        max_homopolymer_length: Maximum junction homopolymer on the full barcode
        spacer: Spacer inserted between the two axis sequences
        junction_role: "left" (column) or "right" (row) junction budget
        partner_right_edge: Right-edge runs of the column codes, for role "right"
        downstream_start_runs: Leading runs of sequences appended after this code
        candidate_mode: "random" or "flip" (see AxisCodeGenerator)
        flip_counts: Bits flipped per "flip" candidate
        require_tolerance: Backtrack instead of keeping an out-of-tolerance code
        
    Returns:
        Tuple of (codes, verification_results)
    """
    generator = AxisCodeGenerator(
        n_positions=n_positions,
        bit_length=bit_length,
        target_distances=target_distances,
        far_features_min_distance=far_features_min_distance,
        max_attempts_per_position=max_attempts_per_position,
        backtrack_depth=backtrack_depth,
        max_backtracks=max_backtracks,
        tolerance=tolerance,
        random_seed=random_seed,
        encoder=encoder,
        g_content_max=g_content_max,
        max_homopolymer_length=max_homopolymer_length,
        spacer=spacer,
        junction_role=junction_role,
        partner_right_edge=partner_right_edge,
        downstream_start_runs=downstream_start_runs,
        candidate_mode=candidate_mode,
        flip_counts=flip_counts,
        require_tolerance=require_tolerance,
    )
    
    codes = generator.generate()
    verification = generator.verify_constraints()
    
    return codes, verification


def generate_column_and_row_codes(
    n_columns: int,
    n_rows: int,
    x_bit_length: int,
    y_bit_length: int,
    target_distances: List[int],
    far_features_min_distance: int,
    encoder,
    g_content_max: float,
    max_homopolymer_length: int,
    spacer: str = "",
    enforce_junction: bool = True,
    max_attempts_per_position: int = 10000,
    backtrack_depth: int = 5,
    max_backtracks: int = 100,
    tolerance: float = 0.5,
    random_seed: Optional[int] = None,
    candidate_mode: str = "random",
    flip_counts: Tuple[int, ...] = (1, 2, 3),
    require_tolerance: bool = False,
    downstream_start_runs: Optional[Dict[str, int]] = None,
) -> Tuple[List[List[int]], List[List[int]], List[str], List[str], Dict, Dict]:
    """
    Generate column codes, then row codes, with composition enforced.

    G content is capped on each axis sequence inside the candidate loop.
    When enforce_junction is true, column codes leave at least one legal
    start base and row codes are rejected when they would extend a column
    edge into a junction homopolymer longer than max_homopolymer_length.
    downstream_start_runs, if set, constrains the trailing run of each row
    code so a later appended sequence (e.g. HSDB bar_2) cannot extend it
    past max_homopolymer_length.
    """
    shared = dict(
        target_distances=target_distances,
        far_features_min_distance=far_features_min_distance,
        max_attempts_per_position=max_attempts_per_position,
        backtrack_depth=backtrack_depth,
        max_backtracks=max_backtracks,
        tolerance=tolerance,
        encoder=encoder,
        g_content_max=g_content_max,
        random_seed=random_seed,
        candidate_mode=candidate_mode,
        flip_counts=flip_counts,
        require_tolerance=require_tolerance,
    )
    left = dict(shared)
    if enforce_junction:
        left.update(
            max_homopolymer_length=max_homopolymer_length,
            spacer=spacer or "",
            junction_role="left",
        )
    x_codes, x_verification = generate_axis_codes(
        n_positions=n_columns,
        bit_length=x_bit_length,
        **left,
    )
    x_sequences = [encoder.encode(code) for code in x_codes]

    right = dict(shared)
    if enforce_junction or downstream_start_runs is not None:
        right["max_homopolymer_length"] = max_homopolymer_length
    if enforce_junction:
        right.update(
            spacer=spacer or "",
            junction_role="right",
            partner_right_edge=max_right_edge_by_base(x_sequences, spacer or ""),
        )
    if downstream_start_runs is not None:
        right["downstream_start_runs"] = downstream_start_runs
    y_codes, y_verification = generate_axis_codes(
        n_positions=n_rows,
        bit_length=y_bit_length,
        **right,
    )
    y_sequences = [encoder.encode(code) for code in y_codes]
    return x_codes, y_codes, x_sequences, y_sequences, x_verification, y_verification
