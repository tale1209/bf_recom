#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
from scipy.stats import binomtest

def binomial_test(k, n, p, alternative, ci_level=0.95):
    """
    Perform exact binomial test
    Parameters:
        k (int): Number of successes
        n (int): Total number of trials
        p (float): Hypothesized success probability
        alternative (str): Test type ('two-sided', 'less', 'greater')
        ci_level (float): Confidence level (default 0.95)
    Returns:
        binomtest result object
    """
    result = binomtest(k, n, p, alternative=alternative)
    return result

def main():
    parser = argparse.ArgumentParser(description="Exact Binomial Test Tool")
    parser.add_argument("-k", "--successes", type=int, required=True, 
                        help="Observed number of successes")
    parser.add_argument("-n", "--trials", type=int, required=True, 
                        help="Total number of trials")
    parser.add_argument("-p", "--prob", type=float, default=0.5, 
                        help="Hypothesized success probability (default 0.5)")
    parser.add_argument("-a", "--alternative", type=str, default="two-sided",
                        choices=["two-sided", "less", "greater"],
                        help="Test alternative: 'two-sided', 'less', 'greater' (default two-sided)")
    parser.add_argument("--ci", type=float, default=0.95,
                        help="Confidence interval level (default 0.95)")

    args = parser.parse_args()
    
    # Validate parameters
    if not 0 <= args.prob <= 1:
        raise ValueError("Hypothesized probability p must be in [0, 1]")
    if args.successes > args.trials or args.successes < 0:
        raise ValueError("Successes cannot exceed trials and must be non-negative")

    result = binomial_test(
        k=args.successes,
        n=args.trials,
        p=args.prob,
        alternative=args.alternative,
        ci_level=args.ci
    )
    
    #Output
    print("===== Binomial Test Results =====")
    print(f"Hypothesis: H0: p = {args.prob}  vs  H1: p {args.alternative}")
    print(f"Observations: Successes = {args.successes}, Trials = {args.trials}")
    print(f"P-value:{result.pvalue:.6f}")
    print(f"Estimated proportion = {result.proportion_estimate:.4f}")
    
    ci_low, ci_high = result.proportion_ci(confidence_level=args.ci)
    print(f"{int(args.ci*100)}% Confidence Interval: [{ci_low:.4f}, {ci_high:.4f}]")
    
    # Statistical conclusion
    alpha = 1 - args.ci
    if result.pvalue < alpha:
        print(f"Conclusion: Reject null hypothesis at α={alpha:.2f} (p < {alpha})")
    else:
        print(f"Conclusion: Fail to reject null hypothesis (p > {alpha})")

if __name__ == "__main__":
    main()