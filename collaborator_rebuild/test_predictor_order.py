"""Check physical label matching rather than histogram axis numbering."""
from collaborator_rebuild.audit_predictor_order import compare_named_atoms


def test_reversed_component_name_order_is_the_same_named_component():
    first = {'atoms': [{'kind': 'R', 'predictors': ['blue', 'core'], 'fraction': .8},
                       {'kind': 'U', 'predictors': ['red'], 'fraction': .2}]}
    reordered = {'atoms': [{'kind': 'U', 'predictors': ['red'], 'fraction': .2},
                           {'kind': 'R', 'predictors': ['core', 'blue'], 'fraction': .8}]}
    assert compare_named_atoms(first, reordered) == (0., 0.)


def test_changed_named_allocation_is_detected():
    first = {'atoms': [{'kind': 'U', 'predictors': ['blue'], 'fraction': .8},
                       {'kind': 'U', 'predictors': ['core'], 'fraction': .2}]}
    changed = {'atoms': [{'kind': 'U', 'predictors': ['core'], 'fraction': .8},
                         {'kind': 'U', 'predictors': ['blue'], 'fraction': .2}]}
    tv, maximum = compare_named_atoms(first, changed)
    assert abs(tv - .6) < 1e-12
    assert abs(maximum - .6) < 1e-12
