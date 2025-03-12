from sympy import *
from sympy.parsing.mathematica import parse_mathematica
from sympy.parsing.latex import parse_latex
from sympy import mathematica_code

import itertools
import numpy as np

def weyl_character(dynkin_labels, fundamental_weights, roots):
    """
    Compute the Weyl character for a given highest weight specified by Dynkin labels.
    
    Parameters:
    - dynkin_labels: List of Dynkin labels [lambda_1, lambda_2, ..., lambda_r]
    - fundamental_weights: List of fundamental weights as numpy arrays
    - roots: List of roots as numpy arrays (used for Weyl group reflections)
    
    Returns:
    - Character function (as a dictionary of weights with coefficients)
    """
    rank = len(dynkin_labels)
    
    # Compute highest weight in terms of fundamental weights
    highest_weight = sum(l * fw for l, fw in zip(dynkin_labels, fundamental_weights))
    
    # Compute Weyl vector (sum of fundamental weights)
    rho = sum(fundamental_weights)
    
    # Weyl group: Generate all sign permutations of coordinate indices
    weyl_group = list(itertools.permutations(range(rank)))
    
    # Compute Weyl alternation sum
    char_dict = {}
    for perm in weyl_group:
        sign = (-1) ** (len(list(itertools.permutations(perm))) % 2)  # Signature of permutation
        w_lambda_rho = np.array([highest_weight[i] + rho[i] for i in perm])
        weight_tuple = tuple(w_lambda_rho - rho)  # Get the weight after shifting back
        char_dict[weight_tuple] = char_dict.get(weight_tuple, 0) + sign
    
    return char_dict

'''# Example for SU(3)
fundamental_weights_su3 = [np.array([2/3, -1/3]), np.array([1/3, 1/3])]
roots_su3 = [np.array([1, -1]), np.array([1, 0]), np.array([0, 1])]  # Example roots

dynkin_labels = [1, 0]  # (1,0) fundamental representation of SU(3)
char = weyl_character(dynkin_labels, fundamental_weights_su3, roots_su3)
print("Character weights:", char)
'''

def PE(polynomial) -> 'polynomial':
    '''
    Plythestic exponetial (PE) of multivariate polynomial
    return: the evaluated form as polynomial
    '''

    # list of monomials
    monomials = expand(polynomial).args

    # factor out the coefficient of each monomial List[Tuple(coefficient, monic monomial)]
    coefficients_monicmonomials = [monomial.as_coeff_Mul() for monomial in monomials]

    # using product-sum formula of PE
    result = 1
    for i in coefficients_monicmonomials:
        term = (1-i[1])**(-i[0])
        result = Mul(result, term)

    return result

def fund_to_anti(character) -> 'polynomial':
    '''
    mapping x_i <-> x_i^{-1} 
    return: the character in antifundamental representation
    '''

    # list of monomials
    monomials = expand(character).args

    # factor out the coefficient of each monomial List[Tuple(coefficient, monic monomial)]
    coefficients_monicmonomials = [monomial.as_coeff_Mul() for monomial in monomials]

    result = 0
    for i in coefficients_monicmonomials:
        term = i[0] * i[1]**(-1)
        result = Add(result, term)
    
    return result



w1, w2, w3, x1, x2, x3, x4, t = symbols('w1 w2 w3 x1 x2 x3 x4 t')
init_printing(use_unicode=True)


penumeratorinput = (w1 + 1/w1 + x1 + 1/x1 + x2 + 1/x2 + x3 + 1/x3 + x4 + 1/x4)*(w2 + 1/w2 + w3 + 1/w3 + 1)*t

pedenominatorinput = parse_mathematica('(w1^2 + 1 + 1/w1^2) t^2 + (w2 + 1/w2 + w3 + 1/w3 + w2 w3 + 1/(w2 w3) + w2/w3 + w3/w2) t^2')


print(mathematica_code(PE(penumeratorinput)))



