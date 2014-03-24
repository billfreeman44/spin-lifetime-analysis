import sympy
import scipy_data_fitting

class Spintronics(scipy_data_fitting.Model):
    """
    Model which sets up common features as explained in
    section 4.3 (standard form) of the notes.

    Should not be used directly, but only subclassed.
    The subclass must define `self.set_f()` which should set `self.expressions['f']`.
    """

    def __init__(self):
        super().__init__()
        self.set_symbols()
        self.set_replacements()
        self.set_replacement_groups()
        self.set_f()
        self.set_expresions()

    def set_symbols(self):
        self.add_symbols(
            'ζ', 'θ', 'ϕ', 'ν', 'r',
            'B', 'ω', 'τ', 'D', 'λ',
            'L', 'W', 'W_F', 'P', 'p', "p'", 'P_σ', 'P_Σ',
            'R_SQ','R_F', 'R_C', 'Ω_F', 'Ω_C',
            'σ_N', 'σ_G', 'ρ_F', 'λ_F', 'A_J', 'd',
            'μ_B', 'ħ')

    def set_replacements(self):
        s = self.symbol
        sqrt = sympy.functions.sqrt

        W, L, D, λ, B, ω, τ, r, μ_B, ħ = self.get_symbols(
            'W', 'L', 'D', 'λ', 'B', 'ω', 'τ', 'r', 'μ_B', 'ħ')

        W_F, R_SQ, R_F, R_C, Ω_F, Ω_C, σ_N, σ_G, ρ_F, λ_F, A_J, d = self.get_symbols(
            'W_F', 'R_SQ', 'R_F', 'R_C', 'Ω_F', 'Ω_C', 'σ_N', 'σ_G', 'ρ_F', 'λ_F', 'A_J', 'd')

        P_σ, P_Σ = self.get_symbols('P_σ', 'P_Σ')

        self.replacements = {
            'ω': ( s('ω'), 2 * μ_B * B / ħ ),
            'ζ': ( s('ζ'), W / λ ),
            'θ': ( s('θ'), λ / r ),
            'ϕ': ( s('ϕ'), L / λ ),
            'ν': ( s('ν'), ω * τ ),
            'λ': ( s('λ'), sqrt(τ * D) ),

            'r': ( s('r'), W * (R_F + R_C) / R_SQ ),
            'R_SQ': ( s('R_SQ'), W / σ_N ),
            'R_F': ( s('R_F'), W * W_F * Ω_F ),
            'R_C': ( s('R_C'), W * W_F * Ω_C ),
            'σ_N': ( s('σ_N'), σ_G / L ),
            'Ω_F': ( s('Ω_F'), ρ_F * λ_F / A_J ),
            'A_J': ( s('A_J'), W * d),

            'p': ( s('p'), (P_σ * R_F  + P_Σ * R_C) / (R_F + R_C) ),
            "p'": ( s("p'"), P_Σ - (P_σ + P_Σ) * R_F / (R_F + R_C) ),
            'p to P': ( s('p'), s('P') ),
            "p' to P": ( s("p'"), s('P') ),

            'θ_zero': ( s('θ'), 0 ),
            'ν_zero': ( s('ν'), 0 ),
        }

    def set_replacement_groups(self):
        self.replacement_groups = {
            'ratios': ('ζ', 'θ', 'ϕ', 'ν'),
            'resistances': ('r', 'R_F', 'R_C', 'R_SQ', 'σ_N'),
            "p, p' to P": ('p to P', "p' to P"),
        }

    def set_expresions(self):
        s = self.symbol
        ex = self.expressions

        W, L, R_SQ, P, ζ = self.get_symbols('W', 'L', 'R_SQ', 'P', 'ζ')
        p, pp, = self.get_symbols('p', "p'")

        ex['λ'] = s('λ')
        ex['r'] = s('r')
        ex['ζ'] = s('ζ')
        ex['θ'] = s('θ')
        ex['ϕ'] = s('ϕ')
        ex['Ω_F'] = s('Ω_F')
        ex['P'] = sympy.functions.sqrt(p * pp)
        ex['R_N'] = s('λ') / (L * W * s('σ_N'))

        ex['f_max'] = self.replace('f', 'ν_zero')
        ex['g'] = ex['f'] / ex['f_max']
        ex['-g'] = - ex['g']
        ex['abs(g)'] = abs(ex['g'])

        ex['nonlocal_resistance'] = R_SQ * (p * pp / ζ) * ex['f']
        ex['nonlocal_resistance_scaled'] = (W * L)**(-1) * ex['nonlocal_resistance']
        ex['nonlocal_resistance_scaled_antiparallel'] = -1 * ex['nonlocal_resistance_scaled']
        ex['nonlocal_resistance_scaled_difference'] = \
            abs(ex['nonlocal_resistance_scaled'] - ex['nonlocal_resistance_scaled_antiparallel'])

class NonZeroField(Spintronics):
    """
    Model for the fixed current nonzero field solution
    as given in section 4.5.1 of the notes.
    """

    def __init__(self):
        self.name = 'nonzero_field_fixed_current'
        super().__init__()

    def set_f(self):
        i = sympy.I
        fn = sympy.functions
        ν, ϕ, θ = self.get_symbols('ν', 'ϕ', 'θ')

        ex_sqrt = fn.sqrt(1 + i * ν)
        ex_denom_1 = 2 * (ex_sqrt  + θ) * fn.exp(ϕ * ex_sqrt)
        ex_denom_2 =  θ**2 * fn.sinh(ϕ * ex_sqrt) * ex_sqrt**(-1)
        self.expressions['f'] = fn.re( (ex_denom_1 + ex_denom_2)**(-1) )

class TransparentContacts(scipy_data_fitting.Model):
    """
    Model given in equation 3 of PhysRevLett.105.167202 for transparent contacts.
    """

    def __init__(self):
        super().__init__()
        self.name = 'transparent_contacts'
        self.add_symbols(
            'u', 'α', 'β', 'A', 'B', 'a', 'L',
            'S_0', 'D', 'τ', 'λ', 'ω',
            'μ_B', 'ħ')

        s = self.symbol
        fn = sympy.functions
        ex = self.expressions
        u, α, β, A, B, a, L, τ, D = self.get_symbols('u', 'α', 'β', 'A', 'B', 'a', 'L', 'τ', 'D')

        self.replacements['a'] = ( s('a'), 2 * s('μ_B') / s('ħ') )
        self.replacements['τ'] = ( s('τ'), α / a )
        self.replacements['D'] = ( s('D'), L**2 / (4 * β * τ) )

        ex['integrand'] = \
            A * fn.sqrt(u)**(-1) * fn.exp(- u - β / u) * fn.cos(B * α * u)

        ex['integrand, B = 0'] = A * fn.sqrt(u)**(-1) * fn.exp(- u - β / u)

        ex['λ'] = fn.sqrt(τ * D)
        ex['L/λ'] = L / ex['λ']
