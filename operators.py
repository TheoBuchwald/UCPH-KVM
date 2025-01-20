from copy import deepcopy

VIR = {"a", "b", "c", "d", "e", "f", "g", "h", "q", "r", "t", "w", "y"}
OCC = {"i", "j", "k", "l", "m", "n", "o", "p", "s", "u", "v", "x", "z"}

class F:
    def __init__(self, indices: list[str], t1_transformed: bool) -> None:
        """Fock operator class."""
        assert len(indices) == 2, f"The number of indices in a Fock operator must be equal to two. Here it was {len(indices)=}."
        self.indices = indices
        self.t1_transformed = t1_transformed
        self.type = "F"

    def __eq__(self, __value: object) -> bool:
        """Check if Fock integral is equal to other operator/amplitude."""
        if type(__value) != type(self):
            return False
        if self.t1_transformed != __value.t1_transformed:
            return False

        p,q = self.indices
        if [p,q] == __value.indices:
            return True
        if (
            [q,p] == __value.indices
            and not self.t1_transformed
        ):
            return True
        return False

    def __str__(self) -> str:
        """String representation of F."""
        p,q = self.indices
        return f"F_{{{p}{q}}} "

    def update_indices(self, old_indices: list[str], new_indices: list[str], translate=False):
        if not translate:
            for o, n in zip(old_indices, new_indices):
                assert o[0] == n[0], f"Old and new index must both be of the same type - here they were {o[0]} and {n[0]}."
        old = []
        new = []
        for o, n in zip(old_indices, new_indices):
            if o in self.indices:
                old.append(o)
                new.append(n)
        replace_indices = []
        for o in old:
            replace_indices.append(self.indices.index(o))
        indices = deepcopy(self.indices)
        for replace, n in zip(replace_indices, new):
            indices[replace] = n
        return F(indices, self.t1_transformed)

    @property
    def dims(self):
        _dims = ""
        for idx in self.indices:
            if idx.lower() in VIR:
                _dims += "v"
            elif idx.lower() in OCC:
                _dims += "o"
            else:
                _dims += idx[0]
        return _dims


class g:
    def __init__(self, indices: list[str], t1_transformed: bool) -> None:
        """Integral class."""
        assert len(indices) == 4, f"The number of indices in an integral must be four. Here it was {len(indices)=}."
        self.indices = indices
        self.t1_transformed = t1_transformed
        self.type = "g"

    def __eq__(self, __value: object) -> bool:
        """Check if integral is equal to other operator/amplitude."""
        if type(__value) != type(self):
            return False
        if self.t1_transformed != __value.t1_transformed:
            return False

        p,q,r,s = self.indices
        if [p,q,r,s] == __value.indices:
            return True
        if [r,s,p,q] == __value.indices:
            return True

        if not self.t1_transformed:
            if ([p,q,s,r] == __value.indices):
                return True
            if ([q,p,r,s] == __value.indices):
                return True
            if ([q,p,s,r] == __value.indices):
                return True
            if ([r,s,q,p] == __value.indices):
                return True
            if ([s,r,p,q] == __value.indices):
                return True
            if([s,r,q,p] == __value.indices):
                return True
        return False

    def __str__(self) -> str:
        """String representation of g."""
        p,q,r,s = self.indices
        return f"g_{{{p}{q}{r}{s}}} "

    def update_indices(self, old_indices: list[str], new_indices: list[str], translate=False):
        if not translate:
            for o, n in zip(old_indices, new_indices):
                assert o[0] == n[0], f"Old and new index must both be of the same type - here they were {o[0]} and {n[0]}."
        old = []
        new = []
        for o, n in zip(old_indices, new_indices):
            if o in self.indices:
                old.append(o)
                new.append(n)
        replace_indices = []
        for o in old:
            replace_indices.append(self.indices.index(o))
        indices = deepcopy(self.indices)
        for replace, n in zip(replace_indices, new):
            indices[replace] = n
        return g(indices, self.t1_transformed)

    @property
    def dims(self):
        _dims = ""
        for idx in self.indices:
            if idx.lower() in VIR:
                _dims += "v"
            elif idx.lower() in OCC:
                _dims += "o"
            else:
                _dims += idx[0]
        return _dims


class L:
    def __init__(self, indices: list[str], t1_transformed: bool) -> None:
        """L integral class."""
        assert len(indices) == 4, f"The number of indices in a L must be equal to four. Here it was {len(indices)=}."
        p,q,r,s = indices
        self.g1 = g([p,q,r,s],t1_transformed)
        self.g2 = g([p,s,r,q],t1_transformed)
        self.type = "L"

    def __eq__(self, __value: object) -> bool:
        """Check if L integral is equal to other operator/amplitude."""
        if type(__value) != type(self):
            return False
        return (
            self.g1 == __value.g1
            and self.g2 == __value.g2
        )

    def __str__(self) -> str:
        """String representation of L."""
        p,q,r,s = self.g1.indices
        return f"L_{{{p}{q}{r}{s}}} "

    def update_indices(self, old_indices: list[str], new_indices: list[str], translate=False):
        if not translate:
            for o, n in zip(old_indices, new_indices):
                assert o[0] == n[0], f"Old and new index must both be of the same type - here they were {o[0]} and {n[0]}."
        old = []
        new = []
        for o, n in zip(old_indices, new_indices):
            if o in self.g1.indices:
                old.append(o)
                new.append(n)
        replace_indices = []
        for o in old:
            replace_indices.append(self.g1.indices.index(o))
        indices = deepcopy(self.g1.indices)
        for replace, n in zip(replace_indices, new):
            indices[replace] = n
        t1_transformed = self.g1.t1_transformed
        return L(indices, t1_transformed)

    @property
    def indices(self):
        return self.g1.indices

    @property
    def dims(self):
        _dims = ""
        for idx in self.indices:
            if idx.lower() in VIR:
                _dims += "v"
            elif idx.lower() in OCC:
                _dims += "o"
            else:
                _dims += idx[0]
        return _dims

