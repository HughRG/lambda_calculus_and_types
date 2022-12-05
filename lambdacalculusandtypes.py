"""A module to aid in Lambda Calculus and Types computations."""

from typing import Union, Optional, Iterable, Mapping
from frozendict import frozendict
from abc import ABC, abstractmethod


class LambdaTerm:
    """A term from Lambda calculus.

    === Public Attributes ===
    variable:
        If self is a variable, then this will be that variable as a string.
        If self is an application, then this will be the empty string.
        If self is an abstraction, then this will be the abstracted variable as a string.
    left_subterm:
        If self is a variable, then this will be None.
        If self is an application, then this will be as expected.
        If self is an abstraction, then this will be the largest proper subterm.
    right_subterm:
        If self is a variable or abstraction, then this will be None.
        If self is an application, then this will be as expected.
    """
    variable: str
    left_subterm: Optional['LambdaTerm']
    right_subterm: Optional['LambdaTerm']

    def __init__(self, variable: str = '',
                 left_subterm: Optional['LambdaTerm'] = None, right_subterm: Optional['LambdaTerm'] = None) -> None:
        """Initialize this LambdaTerm.

        :param variable: If self is to be a variable, then this should be that variable as a string.
        If self is to be an application, then this should be the empty string.
        If self is to be an abstraction, then this should be the abstracted variable as a string.
        :param left_subterm: If self is to be a variable, then this should be None.
        If self is to be an application, then this should be as expected.
        If self is to be an abstraction, then this should be the largest proper subterm.
        :param right_subterm: If self is to be a variable or abstraction, then this should be None.
        If self is to be an application, then this should be as expected.
        """
        # variable or abstraction
        if variable:
            if right_subterm is not None:
                raise ValueError('The inputs do not correspond to a variable, application, or abstraction.')
        # application
        elif left_subterm is None or right_subterm is None:
            raise ValueError('The inputs do not correspond to a variable, application, or abstraction.')

        self.variable = variable
        self.left_subterm = left_subterm
        self.right_subterm = right_subterm

    def __str__(self):
        r"""A string respresentation of this lambda term that can be parsed by LaTeX.

        :return: a string respresentation of this lambda term that can be parsed by LaTeX

        >>> l = LambdaTerm.string_to_lambda_term(r'\lambda y.y')
        >>> print(l)
        \lambda y.y
        >>> l = LambdaTerm.string_to_lambda_term(r'(\lambda y.y)x')
        >>> print(l)
        (\lambda y.y)x
        >>> l = LambdaTerm.string_to_lambda_term(r'\lambda xy.xyy')
        >>> print(l)
        \lambda xy.xyy
        >>> l = LambdaTerm.string_to_lambda_term(r'\lambda xyz.xyy')
        >>> print(l)
        \lambda xyz.xyy
        >>> l = LambdaTerm.string_to_lambda_term(r'x(\lambda x.x)\lambda x.x')
        >>> print(l)
        x(\lambda x.x)(\lambda x.x)
        """
        if self.is_variable():
            return self.variable

        string_representation_of_first_term: str = str(self.left_subterm)

        if self.is_abstraction():
            if self.left_subterm.is_abstraction():
                return fr'\lambda {self.variable}{string_representation_of_first_term[8:]}'
            return fr'\lambda {self.variable}.{self.left_subterm}'

        # application
        if self.left_subterm.is_abstraction():
            string_representation_of_first_term = f'({string_representation_of_first_term})'
        if not self.right_subterm.is_variable():
            return string_representation_of_first_term + f'({self.right_subterm})'
        return string_representation_of_first_term + str(self.right_subterm)

    def __eq__(self, other) -> bool:
        """True iff other is a term that is equal to self under alpha conversion.

        :param other: an object
        :return: True iff other is a term that is equal to self under alpha conversion.

        >>> LambdaTerm.string_to_lambda_term(r'\lambda x.x') == LambdaTerm.string_to_lambda_term(r'\lambda y.y')
        True
        >>> LambdaTerm.string_to_lambda_term(r'\lambda x.x') == LambdaTerm.string_to_lambda_term(r'\lambda y.x')
        False
        """
        if not isinstance(other, LambdaTerm):
            return False

        self_free_variables: frozenset[str] = self.free_variables()
        if other.free_variables() != self_free_variables:
            return False

        if self.is_variable() or other.is_variable():
            return self.is_variable() and other.is_variable()

        if self.is_application() or other.is_application():
            if not self.is_application() and other.is_application():
                return False
            return self.left_subterm == other.left_subterm and self.right_subterm == other.right_subterm

        if self.variable == other.variable:
            return self.left_subterm == other.left_subterm

        other_converted_subterm = other.left_subterm.substitute(LambdaTerm(variable=self.variable), other.variable)

        return self.left_subterm == other_converted_subterm

    @classmethod
    def string_to_lambda_term(cls, lambda_term_string: str) -> 'LambdaTerm':
        r"""Given a string representing a lambda term, return the represented term.

        :param lambda_term_string: A string representing a lambda term where variables are one character
        and there are no extra parentheses.
        :return: The LambdaTerm represented by lambda_term_string.

        >>> l = LambdaTerm.string_to_lambda_term(r'\lambda y.y')
        >>> print(l)
        \lambda y.y
        >>> l = LambdaTerm.string_to_lambda_term(r'(\lambda y.y)x')
        >>> print(l)
        (\lambda y.y)x
        """
        # if it's a variable
        if len(lambda_term_string) == 1:
            return cls(variable=lambda_term_string)

        # if it's an abstraction
        if lambda_term_string.startswith(r'\lambda '):
            end_of_abstracted_variables = lambda_term_string.find('.')
            if end_of_abstracted_variables == 9:
                return cls(variable=lambda_term_string[8],
                           left_subterm=cls.string_to_lambda_term(lambda_term_string[10:]))
            return cls(variable=lambda_term_string[8],
                       left_subterm=cls.string_to_lambda_term(fr'\lambda {lambda_term_string[9:]}'))

        # if it's an application
        application_string_list: list[str] = cls.split_lambda_term_string_into_applications(lambda_term_string)
        application_list: list[cls] = [cls.string_to_lambda_term(term_string)
                                       for term_string in application_string_list]
        application_list.reverse()
        current_term: cls = application_list.pop()
        while application_list:
            current_term = cls(left_subterm=current_term, right_subterm=application_list.pop())
        return current_term

    @classmethod
    def split_lambda_term_string_into_applications(cls, lambda_term_string: str) -> list[str]:
        """Given a string representing a lambda term,
        it maximally splits the string into a list of the form [s, t, u, v, w],
        where each entry is a string representing a subterm of the original term,
        and such that for s', t', u', v', w' being the terms represented by the entries,
        the original term equals s't'u'v'w'

        :param lambda_term_string: a string representing a lambda term
        :return: if lambda_term_string is of the form stuvw it returns [s, t, u, v, w]
        """
        if not lambda_term_string:
            return []

        # if the term is just one abstraction
        if lambda_term_string.startswith(r'\lambda '):
            return [lambda_term_string]

        # if it starts with a variable
        if not lambda_term_string.startswith('('):
            return [lambda_term_string[0]] + cls.split_lambda_term_string_into_applications(lambda_term_string[1:])

        # if it starts with an open parenthesis
        end_of_first_term: int = 0
        open_parentheses: int = 1
        while open_parentheses:
            end_of_first_term += 1
            if lambda_term_string[end_of_first_term] == '(':
                open_parentheses += 1
            elif lambda_term_string[end_of_first_term] == ')':
                open_parentheses -= 1
        return ([lambda_term_string[1:end_of_first_term]]
                + cls.split_lambda_term_string_into_applications(lambda_term_string[end_of_first_term + 1:]))

    def is_variable(self) -> bool:
        """:return: True iff it is a variable."""
        return self.left_subterm is None

    def is_application(self) -> bool:
        """:return: True iff it is an application."""
        return not self.variable

    def is_abstraction(self) -> bool:
        """:return: True iff it is an abstraction."""
        return self.variable and self.left_subterm is not None

    def free_variables(self) -> frozenset[str]:
        """:return: the free variables in the term

        >>> LambdaTerm(variable='x').free_variables()
        frozenset({'x'})
        >>> LambdaTerm.string_to_lambda_term(r'\lambda x.zxzy').free_variables() == frozenset({'y', 'z'})
        True
        >>> LambdaTerm.string_to_lambda_term(r'\lambda n.(\lambda m.(\lambda nmf.n(mf))mm)'
        ...                       + r'((\lambda nfx.n(\lambda th.h(tf))(\lambda g.x)(\lambda y.y))n)').free_variables()
        frozenset()
        """
        if self.is_variable():
            return frozenset(self.variable)
        if self.is_application():
            return self.left_subterm.free_variables() | self.right_subterm.free_variables()
        return self.left_subterm.free_variables() - frozenset(self.variable)

    def is_beta_redex(self) -> bool:
        """:return: True iff if self is a beta redex."""
        return self.is_application() and self.left_subterm.is_abstraction()

    def is_beta_normal_form(self) -> bool:
        """:return: True iff if self is in beta normal form."""
        if self.is_variable():
            return True
        if self.is_abstraction():
            return self.left_subterm.is_beta_normal_form()
        if self.is_beta_redex():
            return False
        # It's an application
        return self.left_subterm.is_beta_normal_form() and self.right_subterm.is_beta_normal_form()

    def substitute(self, term: 'LambdaTerm', variable: str) -> 'LambdaTerm':
        """Return the result of substituting all free occurences of variable with term.

        :param term: to substitute in
        :param variable: to be substituted
        :return: the result of substituting all free occurences of variable with term
        """
        if self.is_variable():
            if self.variable == variable:
                return term
            return self
        if self.is_application():
            return LambdaTerm(left_subterm=self.left_subterm.substitute(term, variable),
                              right_subterm=self.right_subterm.substitute(term, variable))
        if self.variable == variable:
            return self
        return LambdaTerm(variable=self.variable, left_subterm=self.left_subterm.substitute(term, variable))

    def leftmost_beta_reduction(self) -> tuple['LambdaTerm', bool]:
        """:return: a tuple whose first element itself if it's in beta normal form
        and the result of a leftmost reduction otherwise,
        and whose second element is False iff self is in beta normal form

        >>> omega = LambdaTerm.string_to_lambda_term(r'(\lambda x.xx)\lambda x.xx')
        >>> omega_reduct, omega_reduced = omega.leftmost_beta_reduction()
        >>> omega_reduct == omega and omega_reduced
        True
        """
        if self.is_variable():
            return self, False
        if self.is_beta_redex():
            return self.left_subterm.left_subterm.substitute(self.right_subterm, self.left_subterm.variable), True
        reduct, reduced = self.left_subterm.leftmost_beta_reduction()
        if self.is_abstraction():
            return LambdaTerm(variable=self.variable, left_subterm=reduct), reduced
        if reduced:
            return LambdaTerm(left_subterm=reduct, right_subterm=self.right_subterm), True
        reduct, reduced = self.right_subterm.leftmost_beta_reduction()
        return LambdaTerm(left_subterm=self.left_subterm, right_subterm=reduct), reduced

    def is_eta_redex(self) -> bool:
        """:return: True iff if self is a eta redex.

        >>> LambdaTerm.string_to_lambda_term(r'\lambda x.yx').is_eta_redex()
        True
        >>> LambdaTerm.string_to_lambda_term(r'\lambda x.xy').is_eta_redex()
        False
        >>> LambdaTerm.string_to_lambda_term(r'\lambda x.xx').is_eta_redex()
        False
        """
        return (self.is_abstraction() and self.left_subterm.is_application()
                and self.variable == self.left_subterm.right_subterm.variable
                and self.variable not in self.left_subterm.left_subterm.free_variables())

    def is_beta_eta_normal_form(self) -> bool:
        """:return: True iff if self is in beta-eta normal form."""
        if self.is_variable():
            return True
        if self.is_beta_redex() or self.is_eta_redex():
            return False
        return (self.left_subterm.is_beta_eta_normal_form()
                and (self.is_abstraction() or self.right_subterm.is_beta_eta_normal_form()))

    def leftmost_beta_eta_reduction(self) -> tuple['LambdaTerm', bool]:
        """:return: a tuple whose first element itself if it's in beta eta normal form
        and the result of a leftmost beta eta reduction otherwise,
        and whose second element is False iff self is in beta eta normal form

        >>> omega = LambdaTerm.string_to_lambda_term(r'(\lambda x.xx)\lambda x.xx')
        >>> omega_reduct, omega_reduced = omega.leftmost_beta_eta_reduction()
        >>> omega_reduct == omega and omega_reduced
        True
        """
        if self.is_variable():
            return self, False
        if self.is_beta_redex():
            return self.left_subterm.left_subterm.substitute(self.right_subterm, self.left_subterm.variable), True
        if self.is_eta_redex():
            return self.left_subterm.left_subterm, True
        reduct, reduced = self.left_subterm.leftmost_beta_eta_reduction()
        if self.is_abstraction():
            return LambdaTerm(variable=self.variable, left_subterm=reduct), reduced
        if reduced:
            return LambdaTerm(left_subterm=reduct, right_subterm=self.right_subterm), True
        reduct, reduced = self.right_subterm.leftmost_beta_eta_reduction()
        return LambdaTerm(left_subterm=self.left_subterm, right_subterm=reduct), reduced


class ProofTree(ABC):
    """A proof tree that can be printed in LaTeX.

    === Public Attributes ===
    label:
        A label for the last step of this proof tree.
    children:
        The proof trees that start on the line above the base of this proof tree.

    === Representation Invariants ===
    - len(children) <= 2
    """
    label: str
    children: tuple['ProofTree', ...]

    def __init__(self, label: str, children: Iterable['ProofTree']) -> None:
        self.label = label
        self.children = tuple(children)

    def __str__(self) -> str:
        r"""I use the following environment which allows for scaled proof trees,
        as often proof trees are too big to fit on a page.
        The argument of this environment gives the scaled factor (I often use values less than 1).
\newenvironment{scprooftree}[1]
  {\gdef\scalefactor{#1}\begin{center}\proofSkipAmount \leavevmode}
  {\scalebox{\scalefactor}{\DisplayProof}\proofSkipAmount \end{center} }
        """
        return r'\begin{scprooftree}{1}' + '\n' + self._inner_str(1) + '\n' + r'\end{scprooftree}'

    @abstractmethod
    def __hash__(self) -> int:
        pass

    @abstractmethod
    def __eq__(self, other) -> bool:
        pass

    def _inner_str(self, indent: int) -> str:
        """A LaTeX printout of the proof tree at the given indent.

        :param indent: How many indents should precede each line.
        :return: A LaTeX printout of the proof tree at the given indent.
        """
        if not len(self.children) and not self.label:
            return '\t' * indent + r'\AxiomC{\(' + self._root_str() + r'\)}'
        inner_str: str = ''
        for child in self.children:
            inner_str += child._inner_str(indent + 1) + '\n'
        if not len(self.children):
            inner_str += '\t' * (indent + 1) + r'\AxiomC{}' + '\n'
        inner_str += '\t' * indent + r'\LeftLabel{(' + self.label + ')}' + '\n'
        if len(self.children) == 2:
            inner_str += '\t' * indent + r'\BinaryInfC{\('
        else:
            inner_str += '\t' * indent + r'\UnaryInfC{\('
        inner_str += self._root_str() + r'\)}'
        return inner_str

    @abstractmethod
    def _root_str(self) -> str:
        """:return: A string representing the conclusion at the root of the tree."""


class BetaEtaTree(ProofTree):
    """A beta-eta equality proof tree.

    === Public Attributes ===
    l_term:
        The term on the left side of the equality sign at the root.
    r_term:
        The term on the right side of the equality sign at the root.
    """
    l_term: LambdaTerm
    r_term: LambdaTerm
    children: tuple['BetaEtaTree', ...]

    def __init__(self, left_term: LambdaTerm, right_term: LambdaTerm, label: str, children: Iterable['BetaEtaTree']):
        self.left_term = left_term
        self.right_term = right_term
        super().__init__(label, children)

    def __eq__(self, other):
        return (isinstance(other, BetaEtaTree)
                and (self.left_term, self.right_term, self.label, self.children)
                == (other.left_term, other.right_term, other.label, other.children))

    def __hash__(self):
        return hash((self.left_term, self.right_term, self.label, self.children))

    def _root_str(self) -> str:
        return f'{self.left_term} = {self.right_term}'

    @classmethod
    def construct_beta_eta_tree(cls, left_term: LambdaTerm, right_term: LambdaTerm) -> 'BetaEtaTree':
        """Construct a beta eta equality proof tree using user input.

        :param left_term: The term on the left of the equality.
        :param right_term: The term on the right of the equality.
        :return: The beta eta equality proof tree in LaTeX.
        """
        while True:
            rule_number: str = ''
            while rule_number not in ('1', '2', '3', '4', '5', '6', '7', '8'):
                rule_number = input(f'''Which rule will you use to prove the following statement?
{left_term} = {right_term}

1: axiom
2: reflexivity
3: symmetry
4: transitivity
5: application
6: abstraction
7: beta
8: eta
''')
            if rule_number == '1':
                return cls(left_term, right_term, '', tuple())
            elif rule_number == '2':
                if left_term != right_term:
                    print('The two terms are not equal.')
                    continue
                return cls(left_term, right_term, 'refl', tuple())
            elif rule_number == '3':
                return cls(left_term, right_term, 'sym', (cls.construct_beta_eta_tree(right_term, left_term),))
            elif rule_number == '4':
                middle_term = LambdaTerm.string_to_lambda_term(input('What is the middle term?\n'))
                return cls(left_term, right_term, 'trans', (cls.construct_beta_eta_tree(left_term, middle_term),
                                                            cls.construct_beta_eta_tree(middle_term, right_term)))
            elif rule_number == '5':
                if not left_term.is_application() or not right_term.is_application():
                    print('Both terms must be applications.')
                    continue
                return cls(left_term, right_term, 'app',
                           (cls.construct_beta_eta_tree(left_term.left_subterm, right_term.left_subterm),
                            cls.construct_beta_eta_tree(left_term.right_subterm, right_term.right_subterm)))
            elif rule_number == '6':
                if not left_term.is_abstraction() or not right_term.is_abstraction():
                    print('Both terms must be abstractions.')
                    continue
                return cls(left_term, right_term, 'abs', (cls.construct_beta_eta_tree(left_term.left_subterm,
                                                                                      right_term.left_subterm),))
            elif rule_number == '7':
                if not left_term.is_beta_redex():
                    print(f'{left_term} is not a beta redex.')
                    continue
                if left_term.leftmost_beta_reduction() != right_term:
                    print(f'The beta redex {left_term} does not reduce to {right_term}.')
                    continue
                return cls(left_term, right_term, r'\(\beta\)', tuple())
            else:
                if not left_term.is_eta_redex():
                    print(f'{left_term} is not an eta redex.')
                    continue
                if left_term.leftmost_beta_eta_reduction() != right_term:
                    print(f'The eta redex {left_term} does not reduce to {right_term}.')
                    continue
                return cls(left_term, right_term, r'\(\eta\)', tuple())


class CombinatoryLogicTerm:
    """A combinatory logic term.

    === Public Attributes ===
    variable:
        If self is a variable, then this will be that variable as a string.
        If self is the constant K, then this will be the string 'K'.
        If self is the constant S, then this will be the string 'S'.
        If self is an application, then this will be the empty string.
    left_subterm:
        If self is an application, then this will be as expected.
        Otherwise, this will be None.
    right_subterm:
        If self is an application, then this will be as expected.
        Otherwise, this will be None.
    """
    variable: str
    left_subterm: Optional['CombinatoryLogicTerm']
    right_subterm: Optional['CombinatoryLogicTerm']

    def __init__(self, variable: str = '', left_subterm: Optional['CombinatoryLogicTerm'] = None,
                 right_subterm: Optional['CombinatoryLogicTerm'] = None) -> None:
        """Initialize the Combinatory Logic term.

        :param variable: the variable if it is a variable (K and S are not allowed as variables),
        K or S if it's the constant K or S respectively,
        and the empty string if it's an application
        :param left_subterm: the left subterm if it is an application
        :param right_subterm: the right subterm if it is an application

        >>> s = CombinatoryLogicTerm(variable='x')
        >>> print(s)
        x
        """
        if (bool(variable) != (left_subterm is None)) or ((left_subterm is None) != (right_subterm is None)):
            raise ValueError('The inputs do not correspond to a variable, constant, or application.')

        self.variable = variable
        self.left_subterm = left_subterm
        self.right_subterm = right_subterm

    def __str__(self):
        if self.variable:
            return self.variable
        if self.right_subterm.is_application():
            return f'{self.left_subterm}({self.right_subterm})'
        return f'{self.left_subterm}{self.right_subterm}'

    def __hash__(self):
        return hash((self.variable, self.left_subterm, self.right_subterm))

    def __eq__(self, other):
        return (isinstance(other, CombinatoryLogicTerm)
                and ((self.variable, self.left_subterm, self.right_subterm)
                     == (other.variable, other.left_subterm, other.right_subterm)))

    def is_variable(self) -> bool:
        """:return: True iff it is a variable."""
        return bool(self.variable)

    def is_application(self) -> bool:
        """:return: True iff it is an application."""
        return self.left_subterm is not None

    def is_k(self) -> bool:
        """:return: True iff it is the constant K."""
        return self.variable == 'K'

    def is_s(self) -> bool:
        """:return: True iff it is the constant S."""
        return self.variable == 'S'

    def variables(self) -> frozenset[str]:
        """:return: the type variables in the term"""
        if self.is_variable():
            return frozenset(self.variable)
        return self.left_subterm.variables() | self.right_subterm.variables()

    def abstraction(self, variable: str) -> 'CombinatoryLogicTerm':
        """Perform lambda abstraction on a combinatory logic term.

        :param variable: The variable being abstracted
        :return: The combinatory logic term resulting from abstracting variable from this term
        """
        if self.variable == variable:
            return CombinatoryLogicTerm(
                left_subterm=CombinatoryLogicTerm(
                    left_subterm=CombinatoryLogicTerm(variable='S'),
                    right_subterm=CombinatoryLogicTerm(variable='K')), right_subterm=CombinatoryLogicTerm(variable='K'))
        if variable not in self.variables():
            return CombinatoryLogicTerm(left_subterm=CombinatoryLogicTerm(variable='K'), right_subterm=self)
        return CombinatoryLogicTerm(
            left_subterm=CombinatoryLogicTerm(
                left_subterm=CombinatoryLogicTerm(variable='S'),
                right_subterm=self.left_subterm.abstraction(variable)),
            right_subterm=self.right_subterm.abstraction(variable))

    @classmethod
    def string_to_combinatory_logic_term(cls, term_string: str) -> 'CombinatoryLogicTerm':
        """Given a string representing a combinatory logic term, return the represented term.

        :param term_string: A string representing a combinatory logic term with single character variables.
        :return: The represented term.
        """
        if term_string[0] == '(':
            term_string = term_string[1:-1]
        if len(term_string) == 1:
            return cls(variable=term_string)

        start_of_second_subterm: int = -1
        open_parentheses: int = 1 if term_string.endswith(')') else 0
        while open_parentheses:
            start_of_second_subterm -= 1
            if term_string[start_of_second_subterm] == '(':
                open_parentheses -= 1
            elif term_string[start_of_second_subterm] == ')':
                open_parentheses += 1
        return cls(
            left_subterm=cls.string_to_combinatory_logic_term(term_string[:start_of_second_subterm]),
            right_subterm=cls.string_to_combinatory_logic_term(term_string[start_of_second_subterm:]))


class SimpleType:
    """A type in the simple type system.

    === Public Attributes ===
    variable:
        If self is a variable, then this will be that variable as a string.
        Otherwise, this will be the empty string.
    left_subtype:
        If self is an arrow type, then this will be as expected.
        Otherwise, this will be None.
    right_subtype:
        If self is an arrow type, then this will be as expected.
        Otherwise, this will be None.
    """
    variable: str
    left_subtype: Optional['SimpleType']
    right_subtype: Optional['SimpleType']

    def __init__(self, variable: str = '',
                 left_subtype: Optional['SimpleType'] = None, right_subtype: Optional['SimpleType'] = None) -> None:
        """Initialize the simple type.

        :param variable: the variable if it is a variable and the empty string otherwise
        :param left_subtype: the left subterm if it is an arrow type
        :param right_subtype: the right subterm if it is an arrow type

        >>> s = SimpleType(variable='a')
        >>> print(s)
        a
        """
        if (bool(variable) != (left_subtype is None)) or ((left_subtype is None) != (right_subtype is None)):
            raise ValueError('The inputs do not correspond to a variable or arrow type.')

        self.variable = variable
        self.left_subtype = left_subtype
        self.right_subtype = right_subtype

    def __str__(self):
        if self.variable:
            return self.variable
        if self.left_subtype.is_variable():
            return fr'{self.left_subtype}>{self.right_subtype}'
        return fr'({self.left_subtype})>{self.right_subtype}'

    def __hash__(self):
        return hash((self.variable, self.left_subtype, self.right_subtype))

    def __eq__(self, other):
        return isinstance(other, SimpleType) and (self.variable, self.left_subtype, self.right_subtype) == (
            other.variable, other.left_subtype, other.right_subtype)

    def is_variable(self) -> bool:
        """:return: True iff it is a variable."""
        return bool(self.variable)

    def variables(self) -> frozenset[str]:
        """:return: the type variables in the type

        >>> a = SimpleType(variable='a_1')
        >>> a.variables()
        frozenset({'a_1'})
        """
        if self.is_variable():
            return frozenset({self.variable})
        return self.left_subtype.variables() | self.right_subtype.variables()

    @classmethod
    def string_to_type(cls, type_string: str) -> 'SimpleType':
        """Given a string representing a simple type, return the represented type.

        :param type_string: A string representing the type where the only characters are part of variables,
        > meaning an arrow, or parentheses.
        :return: The SimpleType represented by type_string.

        >>> s = SimpleType.string_to_type('a')
        >>> print(s)
        a
        >>> s = SimpleType.string_to_type('(a>a>b)>a>b')
        >>> print(s)
        (a>a>b)>a>b
        >>> s = SimpleType.string_to_type('((a>(a>b))>(a>b))')
        >>> print(s)
        (a>a>b)>a>b
        """
        if '>' not in type_string:
            return cls(variable=type_string)

        topmost_arrow_index: int
        if type_string.startswith('('):
            topmost_arrow_index = 1
            open_parentheses: int = 1
            while open_parentheses:
                if type_string[topmost_arrow_index] == '(':
                    open_parentheses += 1
                elif type_string[topmost_arrow_index] == ')':
                    open_parentheses -= 1
                topmost_arrow_index += 1
            if topmost_arrow_index == len(type_string):
                return cls.string_to_type(type_string[1:-1])
        else:
            topmost_arrow_index = type_string.index('>')

        return cls(left_subtype=cls.string_to_type(type_string[:topmost_arrow_index]),
                   right_subtype=cls.string_to_type(type_string[topmost_arrow_index + 1:]))


class TypeContext:
    r"""A type context.

    === Public Attributes ===
    type_assignments:
        A frozendict which maps variables as strings to simple types as SimpleTypes.

    >>> t = TypeContext({})
    >>> print(t)
    <BLANKLINE>
    >>> t = TypeContext({'x': SimpleType(variable='a')})
    >>> print(t)
    \set{x:a}
    """
    type_assignments: frozendict[str, SimpleType]

    def __init__(self, type_assignments: Mapping[str, SimpleType]) -> None:
        self.type_assignments = frozendict(type_assignments)

    def __eq__(self, other) -> bool:
        return isinstance(other, TypeContext) and self.type_assignments == other.type_assignments

    def __hash__(self) -> int:
        return hash(self.type_assignments)

    def __str__(self) -> str:
        """:return: The type context in LaTeX as it would be displayed in a deduction relation.

        >>> print(TypeContext({'y': SimpleType(variable='b'), 'x': SimpleType(variable='a')}))
        \set{x:a, y:b}
        """
        if not self.type_assignments:
            return ''
        return (r'\set{'
                + ', '.join(f'{subject}:{assignment}' for subject, assignment in sorted(self.type_assignments.items()))
                + '}')

    def __getitem__(self, item) -> SimpleType:
        return self.type_assignments[item]

    def __or__(self, other) -> 'TypeContext':
        if not isinstance(other, TypeContext):
            raise TypeError('A TypeContext can only be unioned with the same type.')
        for variable in self.type_assignments:
            if (variable in other.type_assignments
                    and self.type_assignments[variable] != other.type_assignments[variable]):
                raise ValueError('These TypeContexts are not consistent.')
        return TypeContext(self.type_assignments | other.type_assignments)

    def restriction(self, subjects: Iterable[str]) -> 'TypeContext':
        """Return the restriction of self to subjects.

        :param subjects: the variables to restrict this type context to
        :return: this type context restricted to the given variables
        """
        return TypeContext(
            {subject: self.type_assignments[subject] for subject in subjects if subject in self.type_assignments})

    def set(self, variable: str, simple_type: SimpleType) -> 'TypeContext':
        """Return the result of adding {variable: simple_type} to self.

        :param variable: variable to add to the context
        :param simple_type: type to assign variable
        :return: the new context
        """
        return TypeContext(self.type_assignments.set(variable, simple_type))

    def remove(self, variable: str) -> tuple['TypeContext', SimpleType]:
        """Return the result of removing variable from self and the type variable was assigned.
        raises KeyError if variable is not in self.

        :param variable: variable to remove from the context
        :return: the new context and the type of the variable
        """
        return TypeContext(self.type_assignments.delete(variable)), self.type_assignments[variable]


class DeductionTree(ProofTree):
    """A deduction tree.

    === Public Attributes ===
    type_context:
        The type context at the root.
    lambda_term:
        The lambda term at the root.
    simple_type:
        The deduced type for the lambda term at the root.
    """
    type_context: TypeContext
    lambda_term: LambdaTerm
    simple_type: SimpleType
    children: tuple['DeductionTree', ...]

    def __init__(self, type_context: TypeContext, lambda_term: LambdaTerm, simple_type: SimpleType, label: str,
                 children: Iterable['DeductionTree']):
        self.type_context = type_context
        self.lambda_term = lambda_term
        self.simple_type = simple_type
        super().__init__(label, children)

    def __eq__(self, other):
        return (isinstance(other, DeductionTree)
                and (self.type_context, self.lambda_term, self.simple_type, self.label, self.children)
                == (other.type_context, other.lambda_term, other.simple_type, self.label, other.children))

    def __hash__(self):
        return hash((self.type_context, self.lambda_term, self.simple_type, self.label, self.children))

    def _root_str(self) -> str:
        return fr'{self.type_context} \mapsto {self.lambda_term} : {self.simple_type}'.replace('>', r' \Rightarrow ')

    def type_variables(self) -> frozenset[str]:
        """:return: The type variables occuring in this deduction."""
        type_vars: frozenset[str] = self.simple_type.variables()
        for simple_type in self.type_context.type_assignments.values():
            type_vars |= simple_type.variables()
        for child in self.children:
            type_vars |= child.type_variables()
        return type_vars

    @classmethod
    def principal_type_algorithm(cls, lambda_term: LambdaTerm,
                                 disallowed_type_variables: frozenset[str] = frozenset()) -> Optional['DeductionTree']:
        """
        :param lambda_term:
        :param disallowed_type_variables:
        :return: A principal deduction without any disallowed type variables for lambda_term if it exists

        >>> print(DeductionTree.principal_type_algorithm(LambdaTerm(left_subterm=LambdaTerm(variable='x'),
        ...                                           right_subterm=LambdaTerm(variable='x'))))
        None
        """
        if lambda_term.is_variable():
            next_index: int = 1
            while 'a_{' + str(next_index) + '}' in disallowed_type_variables:
                next_index += 1
            type_variable: SimpleType = SimpleType(variable='a_{' + str(next_index) + '}')
            return cls(TypeContext({lambda_term.variable: type_variable}), lambda_term, type_variable, 'var',
                                 [])

        if lambda_term.is_abstraction():
            sub_deduction: cls = cls.principal_type_algorithm(lambda_term.left_subterm, disallowed_type_variables)
            if sub_deduction is None:
                return None
            if lambda_term.variable in lambda_term.left_subterm.free_variables():
                type_context, abstracted_type = sub_deduction.type_context.remove(lambda_term.variable)
                # return Deduction(type_context, lambda_term,
                #                  SimpleType(left_subtype=abstracted_type, right_subtype=sub_deduction.simple_type),
                #                  'abstraction-main', [sub_deduction])
                return cls(type_context, lambda_term,
                                     SimpleType(left_subtype=abstracted_type, right_subtype=sub_deduction.simple_type),
                                     'abs-main', [sub_deduction])
            next_index: int = 1
            used_type_variables = sub_deduction.type_variables() | disallowed_type_variables
            while 'a_{' + str(next_index) + '}' in used_type_variables:
                next_index += 1
            # return Deduction(sub_deduction.type_context, lambda_term,
            #                  SimpleType(left_subtype=SimpleType(variable=f'a_{next_index}'),
            #                             right_subtype=sub_deduction.simple_type),
            #                  'abstraction-vac', [sub_deduction])
            return cls(sub_deduction.type_context, lambda_term,
                                 SimpleType(left_subtype=SimpleType(variable='a_{' + str(next_index) + '}'),
                                            right_subtype=sub_deduction.simple_type),
                                 'abs-vac', [sub_deduction])

        left_sub_deduction: cls = cls.principal_type_algorithm(lambda_term.left_subterm,
                                                                     disallowed_type_variables)
        if left_sub_deduction is None:
            return None
        right_sub_deduction: cls = cls.principal_type_algorithm(lambda_term.right_subterm,
                                                                      (left_sub_deduction.type_variables()
                                                                       | disallowed_type_variables))
        if right_sub_deduction is None:
            return None

        common_free_variables: list[str] = list(lambda_term.left_subterm.free_variables()
                                                & lambda_term.right_subterm.free_variables())
        left_types: list[SimpleType] = [left_sub_deduction.type_context[var] for var in common_free_variables]
        right_types: list[SimpleType] = [right_sub_deduction.type_context[var] for var in common_free_variables]

        if left_sub_deduction.simple_type.is_variable():
            left_types.append(left_sub_deduction.simple_type)
            used_type_variables = (left_sub_deduction.type_variables()
                                   | right_sub_deduction.type_variables()
                                   | disallowed_type_variables)
            next_index: int = 1
            while 'a_{' + str(next_index) + '}' in used_type_variables:
                next_index += 1
            right_types.append(SimpleType(left_subtype=right_sub_deduction.simple_type,
                                          right_subtype=SimpleType(variable='a_{' + str(next_index) + '}')))
        else:
            left_types.append(left_sub_deduction.simple_type.left_subtype)
            right_types.append(right_sub_deduction.simple_type)

        mgu_domain: frozenset[str] = frozenset()
        for simple_type in left_types:
            mgu_domain |= simple_type.variables()
        for simple_type in right_types:
            mgu_domain |= simple_type.variables()

        unifier: Optional[TypeSubstitution] = TypeSubstitution.sequence_most_general_unifier(left_types, right_types,
                                                                                             ((
                                                                                                          left_sub_deduction.type_variables()
                                                                                                          | right_sub_deduction.type_variables())
                                                                                              - mgu_domain))
        if unifier is None:
            return None

        left_sub_deduction = unifier * left_sub_deduction
        right_sub_deduction = unifier * right_sub_deduction

        return cls(left_sub_deduction.type_context | right_sub_deduction.type_context,
                             LambdaTerm(left_subterm=lambda_term.left_subterm, right_subterm=lambda_term.right_subterm),
                             left_sub_deduction.simple_type.right_subtype, 'app',
                             [left_sub_deduction, right_sub_deduction])


class TypeSubstitution:
    """A generalized type substitution, where we can substitute out any type, not just variables.
    Larger (by containment) types take priority for being substituted out.

    === Public Attributes ===
    type_substitution:
        The type context at the root.
    """
    type_substitution: frozendict[SimpleType, SimpleType]

    def __init__(self, type_substitution: Mapping[Union[str, SimpleType], SimpleType]):
        """Initialize this TypeSubstitution.

        :param type_substitution: A Mapping from strings and SimpleTypes to SimpleTypes.
        If a string is given as a key, it is assumed to represent the corresponding variable type.
        """
        self.type_substitution = frozendict({SimpleType(var): value for var, value in type_substitution.items()
                                             if isinstance(var, str)}
                                            | {simple_type: value for simple_type, value in type_substitution.items()
                                               if isinstance(simple_type, SimpleType)})

    def __mul__(self, other):
        """Apply this substitution to another object.
        Defined for SimpleType, TypeContext, DeductionTree, and TypeSubstitution.

        :param other: The object to apply this substitution to.
        :return: The result of applying self to other.

        >>> U = TypeSubstitution({'a': SimpleType.string_to_type('b>d')})
        >>> A = SimpleType.string_to_type('a>c>c')
        >>> B = SimpleType.string_to_type('(b>d)>a')
        >>> print(U * A)
        (b>d)>c>c
        >>> print(U * B)
        (b>d)>b>d
        """
        if isinstance(other, SimpleType):
            if other in self.type_substitution:
                return self.type_substitution[other]
            if other.is_variable():
                return other
            return SimpleType(left_subtype=self * other.left_subtype, right_subtype=self * other.right_subtype)

        if isinstance(other, TypeContext):
            return TypeContext({var: self * simple_type for var, simple_type in other.type_assignments.items()})

        if isinstance(other, DeductionTree):
            return DeductionTree(self * other.type_context, other.lambda_term, self * other.simple_type, other.label,
                                 (self * child for child in other.children))

        if isinstance(other, TypeSubstitution):
            return TypeSubstitution(
                self.type_substitution
                | {type_var: self * simple_value for type_var, simple_value in other.type_substitution.items()})

        raise TypeError(f'Multiplication is not defined between TypeSubstitution and {type(other)}')

    def restrict(self, variables: frozenset[SimpleType]) -> 'TypeSubstitution':
        """Return self but with the domain restricted to variables.

        :param variables: A superset of the new domain.
        :return: self but with the domain restricted to variables.
        """
        return TypeSubstitution({k: v for k, v in self.type_substitution.items() if k in variables})

    @classmethod
    def most_general_unifier(cls, type_1: SimpleType, type_2: SimpleType) -> Optional['TypeSubstitution']:
        """Return a most general unifier for type_1 and type_2 if one exists,
        with domain and range restricted to the variables occuring in type_1 or type_2.
        If no most general unifier exists, return None.

        :param type_1:
        :param type_2:
        :return: A most general unifier for type_1 and type_2 if one exists,
        with domain and range restricted to the variables occuring in type_1 or type_2.
        If no most general unifier exists, return None.

        >>> A = SimpleType.string_to_type('a>c>c')
        >>> B = SimpleType.string_to_type('(b>d)>a')
        >>> U = TypeSubstitution.most_general_unifier(A, B)
        >>> U * A == U * B
        True
        >>> print(U * A)
        (d>d)>d>d
        >>> A = SimpleType.string_to_type('a>a>a')
        >>> B = SimpleType.string_to_type('b>b')
        >>> U = TypeSubstitution.most_general_unifier(A, B)
        >>> U is None
        True
        """
        if type_1 == type_2:
            return cls({})

        if type_1.is_variable():
            if type_1.variable in type_2.variables():
                return None
            return cls({type_1.variable: type_2})

        if type_2.is_variable():
            if type_2.variable in type_1.variables():
                return None
            return cls({type_2.variable: type_1})

        left_unifier: cls = cls.most_general_unifier(type_1.left_subtype, type_2.left_subtype)
        if left_unifier is None:
            return None
        right_unifier: cls = cls.most_general_unifier(left_unifier * type_1.right_subtype,
                                                      left_unifier * type_2.right_subtype)
        if right_unifier is None:
            return None
        return right_unifier * left_unifier

    @classmethod
    def fresh_most_general_unifier(cls, type_1: SimpleType, type_2: SimpleType,
                                   disallowed_type_variables: frozenset[str]) -> Optional['TypeSubstitution']:
        """Return a most general unifier for type_1 and type_2 if one exists,
        with domain restricted to the variables occuring in type_1 or type_2,
        where the most general unification does not contain the disallowed types.
        If no most general unifier exists, return None.

        :param type_1:
        :param type_2:
        :param disallowed_type_variables: The types that cannot appear in the most general unification.
        :return: A most general unifier for type_1 and type_2 if one exists,
        with domain restricted to the variables occuring in type_1 or type_2,
        where the most general unification does not contain the disallowed types.
        If no most general unifier exists, return None.
        """
        mgu: cls = cls.most_general_unifier(type_1, type_2)
        if mgu is None:
            return None

        current_variables: frozenset[str] = (mgu * type_1).variables()
        problem_variables: frozenset[str] = (mgu * type_1).variables() & disallowed_type_variables
        used_variables: set[str] = set(current_variables | disallowed_type_variables)
        next_index: int = 1
        for type_variable in problem_variables:
            while f'a_{next_index}' in used_variables:
                next_index += 1
            new_variable: str = 'a_{' + str(next_index) + '}'
            mgu = cls({type_variable: SimpleType(variable=new_variable)}) * mgu
            used_variables.add(new_variable)

        return mgu

    @classmethod
    def sequence_most_general_unifier(cls, type_sequence_1: list[SimpleType], type_sequence_2: list[SimpleType],
                                      disallowed_type_variables: frozenset[str] = frozenset()
                                      ) -> Optional['TypeSubstitution']:
        """Return a most general unifier for the sequences type_sequence_1 and type_sequence_2 if one exists,
        with domain restricted to the variables occuring in type_sequence_1 or type_sequence_2,
        where the most general unification does not contain the disallowed types.
        If the sequences are of different length or no most general unifier exists, return None.

        :param type_sequence_1: A sequence of types.
        :param type_sequence_2: A sequence of types.
        :param disallowed_type_variables: The types that cannot appear in the most general unification.
        :return: A most general unifier for the sequences type_sequence_1 and type_sequence_2 if one exists,
        with domain restricted to the variables occuring in type_sequence_1 or type_sequence_2,
        where the most general unification does not contain the disallowed types.
        If the sequences are of different length or no most general unifier exists, return None.
        """
        if len(type_sequence_1) != len(type_sequence_2):
            return None

        concat_type_1 = type_sequence_1[0]
        concat_type_2 = type_sequence_2[0]
        for simple_type in type_sequence_1[1:]:
            concat_type_1 = SimpleType(left_subtype=concat_type_1, right_subtype=simple_type)
        for simple_type in type_sequence_2[1:]:
            concat_type_2 = SimpleType(left_subtype=concat_type_2, right_subtype=simple_type)

        return cls.fresh_most_general_unifier(concat_type_1, concat_type_2, disallowed_type_variables)


def get_term_type_assignment_from_context_or_user(lambda_term: LambdaTerm, type_context: TypeContext) -> SimpleType:
    """Given a lambda term and a type context, deduce the type of the lambda term,
    potentially using input from the user to do so.

    :param type_context:
    :param lambda_term:
    :return: The type of lambda_term.
    """
    if lambda_term.is_variable() and lambda_term.variable in type_context.type_assignments:
        return type_context.type_assignments[lambda_term.variable]
    return SimpleType.string_to_type(input(f'What is the type of {lambda_term}?\n'))


def construct_typing_proof_tree() -> str:
    """

    :return:
    :rtype:
    """
    print('Welcome to the Bussproof Aid.')

    type_context_dict: dict[str, SimpleType] = {}

    variable: str = input('Enter the first variable in the type context or the empty string if the context is empty.\n')
    while variable:
        type_context_dict[variable] = SimpleType.string_to_type(
            input('Enter the type of the variable in the format (a>b)>(c>a)>d\n'))
        variable = input('Enter the next variable in the type context or the empty string if the context is empty.\n')

    starting_type_context: TypeContext = TypeContext(type_context_dict)
    starting_term: LambdaTerm = LambdaTerm.string_to_lambda_term(
        input('Please enter the lambda term you wish to type.\n'))
    starting_type: SimpleType = SimpleType.string_to_type(input('Please enter the type.\n'))

    proof_tree: dict[tuple[TypeContext, LambdaTerm, SimpleType], tuple[
        str, int, list[tuple[TypeContext, LambdaTerm, SimpleType]]]] = {}
    unproven_statements_and_depth: list[tuple[tuple[TypeContext, LambdaTerm, SimpleType], int]] = [
        ((starting_type_context, starting_term, starting_type), 0)]

    while unproven_statements_and_depth:
        statement, depth = unproven_statements_and_depth.pop()
        type_context, lambda_term, simple_type = statement

        if lambda_term.is_variable():
            proof_tree[(type_context, lambda_term, simple_type)] = ('variable', depth, [])
        elif lambda_term.is_application():
            right_subterm_type: SimpleType = get_term_type_assignment_from_context_or_user(lambda_term.right_subterm,
                                                                                           type_context)

            left_subterm_context: TypeContext = type_context.restriction(lambda_term.left_subterm.free_variables())

            right_subterm_context: TypeContext = type_context.restriction(lambda_term.right_subterm.free_variables())

            left_child: tuple[TypeContext, LambdaTerm, SimpleType] = (left_subterm_context, lambda_term.left_subterm,
                                                                      SimpleType(left_subtype=right_subterm_type,
                                                                                 right_subtype=simple_type))

            right_child: tuple[TypeContext, LambdaTerm, SimpleType] = (
                right_subterm_context, lambda_term.right_subterm, right_subterm_type)

            proof_tree[(type_context, lambda_term, simple_type)] = ('application', depth, [right_child, left_child])
            unproven_statements_and_depth.extend([(left_child, depth + 1), (right_child, depth + 1)])

        else:
            assert lambda_term.is_abstraction()
            child_type_context: TypeContext
            abstraction_type: str
            if lambda_term.variable in lambda_term.left_subterm.free_variables():
                child_type_context = type_context.set(lambda_term.variable, simple_type.left_subtype)
                abstraction_type = 'abstraction-main'
            else:
                child_type_context = type_context
                abstraction_type = 'abstraction-vac'
            child: tuple[TypeContext, LambdaTerm, SimpleType] = (
                child_type_context, lambda_term.left_subterm, simple_type.right_subtype)
            proof_tree[(type_context, lambda_term, simple_type)] = (abstraction_type, depth, [child])
            unproven_statements_and_depth.append((child, depth + 1))

    proof_string: str = r'\end{prooftree}'
    unproven_statements: list[tuple[TypeContext, LambdaTerm, SimpleType]] = [
        (starting_type_context, starting_term, starting_type)]
    while unproven_statements:
        statement = unproven_statements.pop()
        type_context, lambda_term, simple_type = statement
        rule, depth, children = proof_tree[statement]
        unproven_statements.extend(children[::-1])

        if rule == 'variable':
            proof_string = ('\t' * depth
                            + r'\UnaryInfC{\(' + fr'{type_context} \mapsto {lambda_term} : {simple_type}' + r'\)}'
                            + f'\n{proof_string}')
            proof_string = '\t' * depth + r'\LeftLabel{(variable)}' + f'\n{proof_string}'
            proof_string = '\t' * (depth + 1) + r'\AxiomC{}' + f'\n{proof_string}'
        elif rule == 'application':
            proof_string = ('\t' * depth
                            + r'\BinaryInfC{\(' + fr'{type_context} \mapsto {lambda_term} : {simple_type}' + r'\)}'
                            + f'\n{proof_string}')
            proof_string = '\t' * depth + r'\LeftLabel{(application)}' + f'\n{proof_string}'
        else:
            proof_string = ('\t' * depth
                            + r'\UnaryInfC{\(' + fr'{type_context} \mapsto {lambda_term} : {simple_type}' + r'\)}'
                            + f'\n{proof_string}')
            proof_string = '\t' * depth + r'\LeftLabel{(' + rule + ')}' + f'\n{proof_string}'
    proof_string = r'\begin{prooftree}' + f'\n{proof_string}'
    proof_string = proof_string.replace('>', r' \Rightarrow ')
    return proof_string


if __name__ == '__main__':
    version = ''
    while version not in ('1', '2'):
        version = input('Enter 1 for beta-eta and 2 for typing.\n')
    if version == '1':
        print('Welcome to the Bussproof Aid')
        l_term: LambdaTerm = LambdaTerm.string_to_lambda_term(input('Please enter the left term.\n'))
        r_term: LambdaTerm = LambdaTerm.string_to_lambda_term(input('Please enter the right term.\n'))
        print(BetaEtaTree.construct_beta_eta_tree(l_term, r_term))
    else:
        print(DeductionTree.principal_type_algorithm(LambdaTerm.string_to_lambda_term(
            input('Please enter the lambda term to type.\n'))))
