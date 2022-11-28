"""A module to aid in Lambda Calculus and Types computations."""

from typing import Union, Optional, Iterable
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

    def __init__(self, variable: str = '', left: Optional['LambdaTerm'] = None, right: Optional['LambdaTerm'] = None):
        """Initialize this LambdaTerm.

        :param variable: If self is to be a variable, then this should be that variable as a string.
        If self is to be an application, then this should be the empty string.
        If self is to be an abstraction, then this should be the abstracted variable as a string.
        :param left: If self is to be a variable, then this should be None.
        If self is to be an application, then this should be as expected.
        If self is to be an abstraction, then this should be the largest proper subterm.
        :param right: If self is to be a variable or abstraction, then this should be None.
        If self is to be an application, then this should be as expected.
        """
        self.variable = variable
        self.left_subterm = left

        # variable or abstraction
        if variable:
            self.right_subterm = None
        # application
        elif left is not None and right is not None:
            self.right_subterm = right
        else:
            raise ValueError('The inputs do not correspond to a variable, application, or abstraction.')

    def __str__(self):
        r"""A string respresentation of this lambda term that can be parsed by LaTeX.

        :return: a string respresentation of this lambda term that can be parsed by LaTeX

        >>> l = string_to_lambda_term(r'\lambda y.y')
        >>> print(l)
        \lambda y.y
        >>> l = string_to_lambda_term(r'(\lambda y.y)x')
        >>> print(l)
        (\lambda y.y)x
        >>> l = string_to_lambda_term(r'\lambda xy.xyy')
        >>> print(l)
        \lambda xy.xyy
        >>> l = string_to_lambda_term(r'\lambda xyz.xyy')
        >>> print(l)
        \lambda xyz.xyy
        >>> l = string_to_lambda_term(r'x(\lambda x.x)\lambda x.x')
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

    def __hash__(self):
        return hash((self.variable, self.left_subterm, self.right_subterm))

    def __eq__(self, other):
        return isinstance(other, LambdaTerm) and (self.variable, self.left_subterm, self.right_subterm) == (
            other.variable, other.left_subterm, other.right_subterm)

    def is_variable(self) -> bool:
        """
        :return: True iff it is a variable.
        """
        return self.left_subterm is None

    def is_application(self) -> bool:
        """
        :return: True iff it is an application.
        """
        return not self.variable

    def is_abstraction(self) -> bool:
        """
        :return: True iff it is an abstraction.
        """
        return self.variable is not None and self.left_subterm is not None

    def free_variables(self) -> frozenset:
        """
        :return: the free variables in the term

        >>> t = string_to_lambda_term(r'\lambda n.(\lambda m.(\lambda nmf.n(mf))mm)((\lambda nfx.n(\lambda th.h(tf))(\lambda g.x)(\lambda y.y))n)')
        >>> t.free_variables()
        frozenset()
        """
        if self.is_variable():
            return frozenset(self.variable)
        if self.is_application():
            return self.left_subterm.free_variables() | self.right_subterm.free_variables()
        return self.left_subterm.free_variables() - frozenset(self.variable)

    def is_beta_redex(self) -> bool:
        """

        :return: True iff if self is a beta redex.
        """
        return self.is_application() and self.left_subterm.is_abstraction()

    def is_beta_normal_form(self) -> bool:
        """

        :return: True iff if self is in beta normal form.
        """
        if self.is_variable():
            return True
        if self.is_abstraction():
            return self.left_subterm.is_beta_normal_form()
        if self.is_beta_redex():
            return False
        # It's an application
        return self.left_subterm.is_beta_normal_form() and self.right_subterm.is_beta_normal_form()

    def substitution(self, term: 'LambdaTerm', variable: str) -> 'LambdaTerm':
        """

        :param term:
        :param variable:
        :return: the result of substituting all free occurences of variable with term
        """
        if variable not in self.free_variables():
            return self
        if self.is_variable():
            return term
        if self.is_application():
            return LambdaTerm(left=self.left_subterm.substitution(term, variable),
                              right=self.right_subterm.substitution(term, variable))
        return LambdaTerm(variable=self.variable, left=self.left_subterm.substitution(term, variable))

    def leftmost_reduction(self) -> 'LambdaTerm':
        """

        :return: self if it's in beta normal form, and the result of a leftmost reduction otherwise.

        >>> omega = string_to_lambda_term(r'(\lambda x.xx)\lambda x.xx')
        >>> omega.leftmost_reduction() == omega
        True
        """
        if self.is_beta_normal_form():
            return self
        if self.is_abstraction():
            return LambdaTerm(variable=self.variable, left=self.left_subterm.leftmost_reduction())
        if self.left_subterm.is_abstraction():
            return self.left_subterm.left_subterm.substitution(self.right_subterm, self.left_subterm.variable)
        if self.left_subterm.is_beta_normal_form():
            return LambdaTerm(left=self.left_subterm, right=self.right_subterm.leftmost_reduction())
        return LambdaTerm(left=self.left_subterm.leftmost_reduction(), right=self.right_subterm)

    def is_eta_redex(self) -> bool:
        """

        :return: True iff if self is a eta redex.
        """
        return (self.is_abstraction() and self.left_subterm.is_application()
                and self.left_subterm.right_subterm.is_variable()
                and self.variable == self.left_subterm.right_subterm.variable)

    def is_beta_eta_normal_form(self) -> bool:
        """

        :return: True iff if self is in beta-eta normal form.
        """
        if self.is_variable():
            return True
        if self.is_beta_redex() or self.is_eta_redex():
            return False
        return (self.left_subterm.is_beta_eta_normal_form()
                and (self.is_abstraction() or self.right_subterm.is_beta_eta_normal_form()))

    def leftmost_beta_eta_reduction(self) -> 'LambdaTerm':
        """

        :return: self if it's in beta eta normal form, and the result of a leftmost beta eta reduction otherwise.
        """
        if self.is_beta_eta_normal_form():
            return self
        if self.is_beta_redex():
            self.left_subterm.substitution(self.right_subterm, self.left_subterm.variable)
        if self.is_eta_redex():
            return self.left_subterm.left_subterm
        if self.is_abstraction():
            return LambdaTerm(variable=self.variable, left=self.left_subterm.leftmost_reduction())
        if self.left_subterm.is_beta_eta_normal_form():
            return LambdaTerm(left=self.left_subterm, right=self.right_subterm.leftmost_reduction())
        return LambdaTerm(left=self.left_subterm.leftmost_reduction(), right=self.right_subterm)


def split_lambda_term_string_into_applications(lambda_term_string: str) -> list[str]:
    """

    :param lambda_term_string: a string representing a lambda term
    :return: if lambda_term_string is of the form stuvw it returns [s, t, u, v, w]
    """
    if not lambda_term_string:
        return []

    # if the term is just one abstraction
    if lambda_term_string.startswith(r'\lambda '):
        return [lambda_term_string]

    # if it starts with a variable
    if lambda_term_string[0] != '(':
        return [lambda_term_string[0]] + split_lambda_term_string_into_applications(lambda_term_string[1:])

    # if it starts with an open parenthesis
    end_of_first_term: int = 0
    open_parentheses: int = 1
    while open_parentheses:
        end_of_first_term += 1
        if lambda_term_string[end_of_first_term] == '(':
            open_parentheses += 1
        elif lambda_term_string[end_of_first_term] == ')':
            open_parentheses -= 1
    return [lambda_term_string[1:end_of_first_term]] + split_lambda_term_string_into_applications(
        lambda_term_string[end_of_first_term + 1:])


def string_to_lambda_term(lambda_term_string: str) -> LambdaTerm:
    r"""

    :param lambda_term_string: A string representing a lambda term where variables are one character
    and there are no extra parentheses.
    :return: The LambdaTerm represented by lambda_term_string.

    >>> l = string_to_lambda_term(r'\lambda y.y')
    >>> print(l)
    \lambda y.y
    >>> l = string_to_lambda_term(r'(\lambda y.y)x')
    >>> print(l)
    (\lambda y.y)x
    """
    # if it's a variable
    if len(lambda_term_string) == 1:
        return LambdaTerm(variable=lambda_term_string)

    # if it's an abstraction
    if lambda_term_string.startswith(r'\lambda '):
        end_of_abstracted_variables = lambda_term_string.find('.')
        if end_of_abstracted_variables == 9:
            return LambdaTerm(variable=lambda_term_string[8],
                              left=string_to_lambda_term(lambda_term_string[10:]))
        return LambdaTerm(variable=lambda_term_string[8],
                          left=string_to_lambda_term(fr'\lambda {lambda_term_string[9:]}'))

    # if it's an application
    application_string_list: list[str] = split_lambda_term_string_into_applications(lambda_term_string)
    application_list: list[LambdaTerm] = [string_to_lambda_term(term_string) for term_string in application_string_list]
    application_list.reverse()
    current_term: LambdaTerm = application_list.pop()
    while application_list:
        current_term = LambdaTerm(left=current_term, right=application_list.pop())
    return current_term


class ProofTree(ABC):
    """
    A proof tree.
    """
    label: str
    children: tuple['ProofTree', ...]

    def __init__(self, label: str, children: Iterable['ProofTree']):
        self.label = label
        self.children = tuple(children)

    def __str__(self):
        return r'\begin{scprooftree}{1}' + '\n' + self._inner_str(1) + '\n' + r'\end{scprooftree}'

    @abstractmethod
    def __hash__(self):
        pass

    @abstractmethod
    def __eq__(self, other):
        pass

    def _inner_str(self, indent: int) -> str:
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
        pass


class BetaEtaTree(ProofTree):
    """
    A beta-eta equality proof tree.
    """
    l_term: LambdaTerm
    r_term: LambdaTerm
    children: tuple['BetaEtaTree', ...]

    def __init__(self, left_term: LambdaTerm, right_term: LambdaTerm, label: str, children: Iterable['BetaEtaTree']):
        self.left_term = left_term
        self.right_term = right_term
        super().__init__(label, children)

    def __hash__(self):
        return hash((self.left_term, self.right_term, self.label, self.children))

    def __eq__(self, other):
        return (isinstance(other, BetaEtaTree)
                and (self.left_term, self.right_term, self.label, self.children)
                == (other.left_term, other.right_term, other.label, other.children))

    def _root_str(self) -> str:
        return f'{self.left_term} = {self.right_term}'


def construct_beta_eta_tree(left_term: LambdaTerm, right_term: LambdaTerm) -> BetaEtaTree:
    """
    Construct a beta eta equality proof tree using user input.

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
            return BetaEtaTree(left_term, right_term, '', tuple())
        elif rule_number == '2':
            if left_term != right_term:
                print('The two terms are not equal.')
                continue
            return BetaEtaTree(left_term, right_term, 'refl', tuple())
        elif rule_number == '3':
            return BetaEtaTree(left_term, right_term, 'sym', (construct_beta_eta_tree(right_term, left_term),))
        elif rule_number == '4':
            middle_term = string_to_lambda_term(input('What is the middle term?\n'))
            return BetaEtaTree(left_term, right_term, 'trans', (construct_beta_eta_tree(left_term, middle_term),
                                                                construct_beta_eta_tree(middle_term, right_term)))
        elif rule_number == '5':
            if not left_term.is_application() or not right_term.is_application():
                print('Both terms must be applications.')
                continue
            return BetaEtaTree(left_term, right_term, 'app',
                               (construct_beta_eta_tree(left_term.left_subterm, right_term.left_subterm),
                                construct_beta_eta_tree(left_term.right_subterm, right_term.right_subterm)))
        elif rule_number == '6':
            if not left_term.is_abstraction() or not right_term.is_abstraction():
                print('Both terms must be abstractions.')
                continue
            return BetaEtaTree(left_term, right_term, 'abs',
                               (construct_beta_eta_tree(left_term.left_subterm, right_term.left_subterm),))
        elif rule_number == '7':
            if not left_term.is_beta_redex():
                print(f'{left_term} is not a beta redex.')
                continue
            if left_term.leftmost_reduction() != right_term:
                print(f'The beta redex {left_term} does not reduce to {right_term}.')
                continue
            return BetaEtaTree(left_term, right_term, r'\(\beta\)', tuple())
        else:
            if not left_term.is_eta_redex():
                print(f'{left_term} is not an eta redex.')
                continue
            if left_term.leftmost_beta_eta_reduction() != right_term:
                print(f'The eta redex {left_term} does not reduce to {right_term}.')
                continue
            return BetaEtaTree(left_term, right_term, r'\(\eta\)', tuple())


class CombinatoryLogicTerm:
    """
    A combinatory logic term.
    """
    variable: Optional[str]
    left_subterm: Optional['CombinatoryLogicTerm']
    right_subterm: Optional['CombinatoryLogicTerm']

    def __init__(self, variable: Optional[str] = None, left_subterm: Optional['CombinatoryLogicTerm'] = None,
                 right_subterm: Optional['CombinatoryLogicTerm'] = None):
        """
        Init should be given the info for exactly one kind of term

        :param variable: the variable if it is a variable, K if it's the constant K, and S if it's the constant S
        :param left_subterm: the left subterm if it is an arrow type
        :param right_subterm: the right subterm if it is an arrow type

        >>> s = CombinatoryLogicTerm(variable='A')
        >>> print(s)
        A
        """
        assert (sum((variable is not None, left_subterm is not None))
                * ((left_subterm is None) == (right_subterm is None)))

        self.variable = variable
        self.left_subterm = left_subterm
        self.right_subterm = right_subterm

    def __str__(self):
        if self.variable is not None:
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
        """
        :return: True iff it is a variable.
        """
        return self.variable is not None

    def is_application(self) -> bool:
        """
        :return: True iff it is an application.
        """
        return self.left_subterm is not None

    def is_k(self) -> bool:
        """
        :return: True iff it is the constant K.
        """
        return self.variable == 'K'

    def is_s(self) -> bool:
        """
        :return: True iff it is the constant S.
        """
        return self.variable == 'S'

    def variables(self) -> frozenset:
        """
        :return: the type variables in the term
        """
        if self.is_variable():
            return frozenset(self.variable)
        return self.left_subterm.variables() | self.right_subterm.variables()

    def combinatory_logic_lambda(self, variable: str) -> 'CombinatoryLogicTerm':
        """
        Perform lambda abstraction on a combinatory logic term.

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
                right_subterm=self.left_subterm.combinatory_logic_lambda(variable)),
            right_subterm=self.right_subterm.combinatory_logic_lambda(variable))


def string_to_combinatory_logic_term(term_string: str) -> CombinatoryLogicTerm:
    """

    :param term_string: A combinatory logic term with single character variables and no extra parentheses
    (other than optionally outer parentheses).
    :return: The combinatory logic term represented by term.
    """
    if term_string[0] == '(':
        term_string = term_string[1:-1]
    if len(term_string) == 1:
        return CombinatoryLogicTerm(variable=term_string)

    start_of_second_subterm: int = -1
    open_parentheses: int = 1 if term_string.endswith(')') else 0
    while open_parentheses:
        start_of_second_subterm -= 1
        if term_string[start_of_second_subterm] == '(':
            open_parentheses -= 1
        elif term_string[start_of_second_subterm] == ')':
            open_parentheses += 1
    return CombinatoryLogicTerm(left_subterm=string_to_combinatory_logic_term(term_string[:start_of_second_subterm]),
                                right_subterm=string_to_combinatory_logic_term(term_string[start_of_second_subterm:]))


class SimpleType:
    """
    A type in the simple type system.
    """
    variable: Optional[str]
    left_subtype: Optional['SimpleType']
    right_subtype: Optional['SimpleType']

    def __init__(self, variable: Optional[str] = None, left: Optional['SimpleType'] = None,
                 right: Optional['SimpleType'] = None):
        """
        Init should be given a variable string XOR a left and right subtype.

        :param variable: the variable if it is a variable
        :param left: the left subtype if it is an arrow type
        :param right: the right subtype if it is an arrow type

        >>> s = SimpleType(variable='A')
        >>> print(s)
        A
        """
        assert (variable is not None or (left is not None and right is not None)
                and (variable is None or (left is None and right is None)))
        if variable is not None:
            self.variable = variable
            self.left_subtype = None
            self.right_subtype = None
        elif left is not None and right is not None:
            self.variable = None
            self.left_subtype = left
            self.right_subtype = right
        else:
            print('Fuck you')
            self.variable = 'A'
            self.left_subtype = None
            self.right_subtype = None

    def __str__(self):
        if self.variable is not None:
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
        """
        :return: True iff it is a variable.
        """
        return self.variable is not None

    def variables(self) -> frozenset[str]:
        """
        :return: the type variables in the type

        >>> a = SimpleType(variable='a_1')
        >>> a.variables()
        frozenset({'a_1'})
        """
        if self.is_variable():
            return frozenset({self.variable})
        return self.left_subtype.variables() | self.right_subtype.variables()


def string_to_type(type_string: str) -> SimpleType:
    """

    :param type_string: A string representing the type where the only characters are part of variables,
    > meaning an arrow, or parentheses, with no extra parentheses (other than optionally outer parentheses).
    :return: The SimpleType represented by type_string.

    >>> s = string_to_type('(A>A>B)>A>B')
    >>> print(s)
    (A>A>B)>A>B
    """
    if type_string.endswith(')'):
        type_string = type_string[1:-1]

    if '>' not in type_string:
        return SimpleType(variable=type_string)

    start_of_second_subtype: int = 1
    open_parentheses: int = 1 if type_string[0] == '(' else 0
    while open_parentheses:
        if type_string[start_of_second_subtype] == '(':
            open_parentheses += 1
        elif type_string[start_of_second_subtype] == ')':
            open_parentheses -= 1
        start_of_second_subtype += 1
    return SimpleType(left=string_to_type(type_string[:start_of_second_subtype]),
                      right=string_to_type(type_string[start_of_second_subtype + 1:]))


class TypeContext:
    r"""
    At least one of indices and extra_assignment should be None.

    Consistency is not checked.

    >>> t = TypeContext({})
    >>> print(t)
    <BLANKLINE>
    >>> t = TypeContext({'x': SimpleType(variable='A')})
    >>> print(t)
    \set{x:A}
    """
    type_assignments: frozendict[str, SimpleType]

    def __init__(self, type_assignments: Union[frozendict[str, SimpleType], dict[str, SimpleType]]):
        if isinstance(type_assignments, frozendict):
            self.type_assignments = type_assignments
        else:
            self.type_assignments = frozendict(type_assignments)

    def __hash__(self):
        return hash(self.type_assignments)

    def __eq__(self, other):
        return isinstance(other, TypeContext) and self.type_assignments == other.type_assignments

    def __str__(self):
        if not self.type_assignments:
            return ''
        subjects = list(self.type_assignments.keys())
        subjects.sort()
        string_repesentation: str = r'\set{'
        for subject in subjects[:-1]:
            string_repesentation += f'{subject}:{self.type_assignments[subject]}, '
        return string_repesentation + f'{subjects[-1]}:{self.type_assignments[subjects[-1]]}' + r'}'

    def __getitem__(self, item):
        return self.type_assignments[item]

    def __or__(self, other):
        if isinstance(other, TypeContext):
            return TypeContext(self.type_assignments | other.type_assignments)
        raise TypeError('A TypeContext can only be unioned with the same type.')

    def restriction(self, subjects: Iterable[str]) -> 'TypeContext':
        """
        :param subjects: the variables to restrict this type context to
        :return: this type context restricted to the given variables
        """
        return TypeContext(
            {subject: self.type_assignments[subject] for subject in subjects if subject in self.type_assignments})

    def set(self, variable: str, simple_type: SimpleType) -> 'TypeContext':
        """
        :param variable: variable to add to the context
        :param simple_type: type to assign variable
        :return: the new context
        """
        return TypeContext(self.type_assignments.set(variable, simple_type))

    def remove(self, variable: str) -> tuple['TypeContext', SimpleType]:
        """
        :param variable: variable to remove from the context
        :return: the type of the variable and the new context
        """
        return TypeContext(self.type_assignments.delete(variable)), self.type_assignments[variable]


class Deduction(ProofTree):
    """
    A typing deduction.
    """
    type_context: TypeContext
    lambda_term: LambdaTerm
    simple_type: SimpleType
    children: tuple['Deduction', ...]

    def __init__(self, type_context: TypeContext, lambda_term: LambdaTerm, simple_type: SimpleType, label: str,
                 children: Iterable['Deduction']):
        self.type_context = type_context
        self.lambda_term = lambda_term
        self.simple_type = simple_type
        super().__init__(label, children)

    def __hash__(self):
        return hash((self.type_context, self.lambda_term, self.simple_type, self.children))

    def __eq__(self, other):
        return (isinstance(other, Deduction)
                and (self.type_context, self.lambda_term, self.simple_type, self.children)
                == (other.type_context, other.lambda_term, other.simple_type, other.children))

    def _root_str(self) -> str:
        return fr'{self.type_context} \mapsto {self.lambda_term} : {self.simple_type}'.replace('>', r' \Rightarrow ')

    def type_variables(self) -> frozenset[str]:
        """
        :return: The type variables occuring in this deduction.
        """
        type_vars: frozenset[str] = self.simple_type.variables()
        for child in self.children:
            type_vars |= child.type_variables()
        return type_vars


class TypeSubstitution:
    """
    A type substitution.
    """
    type_substitution: frozendict[SimpleType, SimpleType]

    def __init__(self, type_substitution: Union[dict[Union[str, SimpleType], SimpleType],
                                                frozendict[Union[str, SimpleType], SimpleType]]):
        """

        :param type_substitution: Must be consistent
        """
        self.type_substitution = frozendict({SimpleType(var): value for var, value in type_substitution.items()
                                             if isinstance(var, str)}
                                            | {simple_type: value for simple_type, value in type_substitution.items()
                                               if isinstance(simple_type, SimpleType)})

    def __mul__(self, other):
        """

        :param other:
        :return:

        >>> U = TypeSubstitution({'a': string_to_type('b>d')})
        >>> A = string_to_type('a>c>c')
        >>> B = string_to_type('(b>d)>a')
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
            return SimpleType(left=self * other.left_subtype, right=self * other.right_subtype)

        if isinstance(other, TypeContext):
            return TypeContext({var: self * simple_type for var, simple_type in other.type_assignments.items()})

        if isinstance(other, Deduction):
            return Deduction(self * other.type_context, other.lambda_term, self * other.simple_type, other.label,
                             (self * child for child in other.children))

        if isinstance(other, TypeSubstitution):
            return TypeSubstitution(
                {type_var: self * simple_value for type_var, simple_value in other.type_substitution.items()}
                | self.type_substitution)

        raise TypeError(f'Multiplication is not defined between TypeSubstitution and {type(other)}')


def get_term_type_assignment_from_context_or_user(type_context: TypeContext, lambda_term: LambdaTerm) -> SimpleType:
    """

    :param type_context:
    :type type_context:
    :param lambda_term:
    :type lambda_term:
    :return:
    :rtype:
    """
    if lambda_term.is_variable() and lambda_term.variable in type_context.type_assignments:
        return type_context.type_assignments[lambda_term.variable]
    return string_to_type(input(f'What is the type of {lambda_term}?\n'))


def construct_typing_proof_tree() -> str:
    """

    :return:
    :rtype:
    """
    print('Welcome to the Bussproof Aid.')

    type_context_dict: dict[str, SimpleType] = {}

    variable: str = input('Enter the first variable in the type context or the empty string if the context is empty.\n')
    while variable:
        type_context_dict[variable] = string_to_type(
            input('Enter the type of the variable in the format (A>B)>(C>A)>D\n'))
        variable = input('Enter the next variable in the type context or the empty string if the context is empty.\n')

    starting_type_context: TypeContext = TypeContext(type_context_dict)
    starting_term: LambdaTerm = string_to_lambda_term(input('Please enter the lambda term you wish to type.\n'))
    starting_type: SimpleType = string_to_type(input('Please enter the type.\n'))

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
            right_subterm_type: SimpleType = get_term_type_assignment_from_context_or_user(type_context,
                                                                                           lambda_term.right_subterm)

            left_subterm_context: TypeContext = type_context.restriction(lambda_term.left_subterm.free_variables())

            right_subterm_context: TypeContext = type_context.restriction(lambda_term.right_subterm.free_variables())

            left_child: tuple[TypeContext, LambdaTerm, SimpleType] = (left_subterm_context, lambda_term.left_subterm,
                                                                      SimpleType(left=right_subterm_type,
                                                                                 right=simple_type))

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


def most_general_unifier(type_1: SimpleType, type_2: SimpleType) -> Optional[TypeSubstitution]:
    """

    :param type_1:
    :param type_2:
    :return: The most general unifier for a and b.

    >>> A = SimpleType(left=SimpleType(variable='a'),
    ...                right=SimpleType(left=SimpleType(variable='c'),
    ...                                         right=SimpleType(variable='c')))
    >>> B = SimpleType(left=SimpleType(left=SimpleType(variable='b'),
    ...                                        right=SimpleType(variable='d')),
    ...                right=SimpleType(variable='a'))
    >>> U = most_general_unifier(A, B)
    >>> U * A == U * B
    True
    >>> print(U * A)
    (d>d)>d>d
    """
    if type_1 == type_2:
        return TypeSubstitution({})

    if type_1.is_variable():
        if type_1.variable in type_2.variables():
            return None
        return TypeSubstitution({type_1.variable: type_2})

    if type_2.is_variable():
        if type_2.variable in type_1.variables():
            return None
        return TypeSubstitution({type_2.variable: type_1})

    left_unifier: TypeSubstitution = most_general_unifier(type_1.left_subtype, type_2.left_subtype)
    if left_unifier is None:
        return None
    right_unifier: TypeSubstitution = most_general_unifier(left_unifier * type_1.right_subtype,
                                                           left_unifier * type_2.right_subtype)
    if right_unifier is None:
        return None
    return right_unifier * left_unifier


def fresh_most_general_unifier(type_1: SimpleType, type_2: SimpleType,
                               disallowed_type_variables: frozenset[str]) -> Optional[TypeSubstitution]:
    """

    :param type_1:
    :param type_2:
    :param disallowed_type_variables:
    :return: an mgu for type 1 and 2 where the images don't intersect the disallowed types
    """
    mgu: TypeSubstitution = most_general_unifier(type_1, type_2)
    if mgu is None:
        return None
    current_variables: frozenset[str] = (mgu * type_1).variables() | (mgu * type_2).variables()
    current_problem_variables: frozenset[str] = current_variables & disallowed_type_variables
    next_index: int = 1
    while current_problem_variables:
        for a in current_problem_variables:
            used_variables: frozenset[str] = current_variables | disallowed_type_variables
            while f'a_{next_index}' in used_variables:
                next_index += 1
            mgu = TypeSubstitution({a: SimpleType(variable='a_{' + str(next_index) + '}')}) * mgu
            current_variables = (mgu * type_1).variables() | (mgu * type_2).variables()
        current_problem_variables = current_variables & disallowed_type_variables

    return mgu


def sequence_most_general_unifier(type_sequence_1: list[SimpleType], type_sequence_2: list[SimpleType],
                                  disallowed_type_variables: frozenset[str] = frozenset()
                                  ) -> Optional[TypeSubstitution]:
    """

    :param type_sequence_1: A sequence of types.
    :param type_sequence_2: A sequence of types.
    :param disallowed_type_variables: Variables that cannot be in the image of the mgu applied to the inputs.
    :return: A most general unifier for the two inputs.
    """
    if len(type_sequence_1) != len(type_sequence_2):
        return None

    used_type_variables: set = set()
    for simple_type in type_sequence_1:
        used_type_variables |= simple_type.variables()
    for simple_type in type_sequence_2:
        used_type_variables |= simple_type.variables()
    next_index: int = -1
    while 'a_{' + str(next_index) + '}' in used_type_variables:
        next_index += 1

    concat_type_1 = SimpleType(variable='a_{' + str(next_index) + '}')
    concat_type_2 = concat_type_1
    for simple_type in type_sequence_1[::-1]:
        concat_type_1 = SimpleType(left=simple_type, right=concat_type_1)
    for simple_type in type_sequence_2[::-1]:
        concat_type_2 = SimpleType(left=simple_type, right=concat_type_2)

    return fresh_most_general_unifier(concat_type_1, concat_type_2, disallowed_type_variables)


def principal_type_algorithm(lambda_term: LambdaTerm,
                             disallowed_type_variables: frozenset[str] = frozenset()) -> Optional[Deduction]:
    """
    :param lambda_term:
    :param disallowed_type_variables:
    :return: A principal deduction without any disallowed type variables for lambda_term if it exists

    >>> print(principal_type_algorithm(LambdaTerm(left=LambdaTerm(variable='x'),
    ...                                           right=LambdaTerm(variable='x'))))
    None
    """
    if lambda_term.is_variable():
        next_index: int = 1
        while 'a_{' + str(next_index) + '}' in disallowed_type_variables:
            next_index += 1
        type_variable: SimpleType = SimpleType(variable='a_{' + str(next_index) + '}')
        # return Deduction(TypeContext({lambda_term.variable: type_variable}), lambda_term, type_variable, 'variable', [])
        return Deduction(TypeContext({lambda_term.variable: type_variable}), lambda_term, type_variable, 'var', [])

    if lambda_term.is_abstraction():
        sub_deduction: Deduction = principal_type_algorithm(lambda_term.left_subterm, disallowed_type_variables)
        if sub_deduction is None:
            return None
        if lambda_term.variable in lambda_term.left_subterm.free_variables():
            type_context, abstracted_type = sub_deduction.type_context.remove(lambda_term.variable)
            # return Deduction(type_context, lambda_term,
            #                  SimpleType(left_subtype=abstracted_type, right_subtype=sub_deduction.simple_type),
            #                  'abstraction-main', [sub_deduction])
            return Deduction(type_context, lambda_term,
                             SimpleType(left=abstracted_type, right=sub_deduction.simple_type),
                             'abs-main', [sub_deduction])
        next_index: int = 1
        used_type_variables = sub_deduction.type_variables() | disallowed_type_variables
        while 'a_{' + str(next_index) + '}' in used_type_variables:
            next_index += 1
        # return Deduction(sub_deduction.type_context, lambda_term,
        #                  SimpleType(left_subtype=SimpleType(variable=f'a_{next_index}'),
        #                             right_subtype=sub_deduction.simple_type),
        #                  'abstraction-vac', [sub_deduction])
        return Deduction(sub_deduction.type_context, lambda_term,
                         SimpleType(left=SimpleType(variable='a_{' + str(next_index) + '}'),
                                    right=sub_deduction.simple_type),
                         'abs-vac', [sub_deduction])

    left_sub_deduction: Deduction = principal_type_algorithm(lambda_term.left_subterm, disallowed_type_variables)
    if left_sub_deduction is None:
        return None
    right_sub_deduction: Deduction = principal_type_algorithm(lambda_term.right_subterm,
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
        right_types.append(SimpleType(left=right_sub_deduction.simple_type,
                                      right=SimpleType(variable='a_{' + str(next_index) + '}')))
    else:
        left_types.append(left_sub_deduction.simple_type.left_subtype)
        right_types.append(right_sub_deduction.simple_type)

    mgu_domain: frozenset[str] = frozenset()
    for simple_type in left_types:
        mgu_domain |= simple_type.variables()
    for simple_type in right_types:
        mgu_domain |= simple_type.variables()

    unifier: Optional[TypeSubstitution] = sequence_most_general_unifier(left_types, right_types,
                                                                        ((left_sub_deduction.type_variables()
                                                                          | right_sub_deduction.type_variables())
                                                                         - mgu_domain))
    if unifier is None:
        return None

    left_sub_deduction = unifier * left_sub_deduction
    right_sub_deduction = unifier * right_sub_deduction

    # return Deduction(left_sub_deduction.type_context | right_sub_deduction.type_context,
    #                  LambdaTerm(left_subterm=lambda_term.left_subterm, right_subterm=lambda_term.right_subterm),
    #                  left_sub_deduction.simple_type.right_subtype, 'application',
    #                  [left_sub_deduction, right_sub_deduction])
    return Deduction(left_sub_deduction.type_context | right_sub_deduction.type_context,
                     LambdaTerm(left=lambda_term.left_subterm, right=lambda_term.right_subterm),
                     left_sub_deduction.simple_type.right_subtype, 'app',
                     [left_sub_deduction, right_sub_deduction])


if __name__ == '__main__':
    version = ''
    while version not in ('1', '2'):
        version = input('Enter 1 for beta-eta and 2 for typing.\n')
    if version == '1':
        print('Welcome to the Bussproof Aid')
        l_term: LambdaTerm = string_to_lambda_term(input('Please enter the left term.\n'))
        r_term: LambdaTerm = string_to_lambda_term(input('Please enter the right term.\n'))
        print(construct_beta_eta_tree(l_term, r_term))
    else:
        print(principal_type_algorithm(string_to_lambda_term(input('Please enter the lambda term to type.\n'))))
