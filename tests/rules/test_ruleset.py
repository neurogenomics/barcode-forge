from barcode_forge.config import Config
from barcode_forge.rules.ruleset import Rule, RuleSet

calls: list[str] = []


def _make(name, cost, result):
    def check(seq, cfg):
        calls.append(name)
        return result

    return Rule(name=name, cost=cost, check=check)


def test_passes_when_all_rules_pass():
    rs = RuleSet([_make("a", 1, True), _make("b", 2, True)], Config())
    assert rs.passes("ACGTACGTAC")
    assert not rs.rejections


def test_runs_cheapest_first_and_stops_at_first_failure():
    calls.clear()
    rs = RuleSet([_make("slow", 10, False), _make("fast", 1, False)], Config())
    assert not rs.passes("ACGTACGTAC")
    assert calls == ["fast"]
    assert rs.rejections == {"fast": 1}


def test_rejections_accumulate():
    rs = RuleSet([_make("a", 1, True), _make("b", 2, False)], Config())
    rs.passes("AAAA")
    rs.passes("CCCC")
    assert rs.rejections == {"b": 2}


def test_names_in_cost_order():
    rs = RuleSet([_make("b", 2, True), _make("a", 1, True)], Config())
    assert rs.names == ("a", "b")
