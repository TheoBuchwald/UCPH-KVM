from collections.abc import Generator
from itertools import permutations
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


class amplitude:
    def __init__(self, indices: list[str]) -> None:
        """Amplitude class."""
        assert len(indices) % 2 == 0, f"The number of indices in an amplitude must be even. Here it was {len(indices)=}."
        self.indices = indices

    def __eq__(self, __value: object) -> bool:
        """Check if amplitude is equal to other amplitude/integral."""
        if type(__value) != type(self):
            return False
        if len(self.indices) != len(__value.indices):
            return False

        vir_index_permutations = permutations(self.indices[::2])
        occ_index_permutations = permutations(self.indices[1::2])

        for vir_idx, occ_idx in zip(vir_index_permutations, occ_index_permutations):
            if (
                vir_idx == tuple(__value.indices[::2])
                and occ_idx == tuple(__value.indices[1::2])
            ):
                return True
        return False

    def __len__(self) -> int:
        return len(self.indices) // 2

    def __str__(self) -> str:
        """String representation of amplitude."""
        string = "t_{"
        for p in self.indices:
            string += p
        string += "} "
        return string

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
        return amplitude(indices)


class t:
    def __init__(self, indices: list[list[str]]) -> None:
        """t class."""
        self.amplitudes: list[amplitude] = []
        for idx in indices:
            self.amplitudes.append(
                amplitude(idx)
            )

    def __str__(self) -> str:
        """String representation of t."""
        string = ""
        for amp in self.amplitudes:
            string += f"{str(amp)}"
        return string

    def __eq__(self, __value: object) -> bool:
        """Check if T is equal to other T."""
        if type(__value) != type(self):
            return False
        for amp1, amp2 in zip(self.amplitudes, __value.amplitudes):
            if amp1 != amp2:
                return False
        return True

    @property
    def indices(self):
        _indices = []
        for amp in self.amplitudes:
            _indices += amp.indices
        return _indices

    def update_indices(self, old_indices: list[str], new_indices: list[str], translate=False):
        if not translate:
            for old, new in zip(old_indices, new_indices):
                assert old[0] == new[0], f"Old and new index must both be of the same type - here they were {old[0]} and {new[0]}."
        t_indices = []
        for amp in self.amplitudes:
            old = []
            new = []
            for o, n in zip(old_indices, new_indices):
                if o in amp.indices:
                    old.append(o)
                    new.append(n)
            replace_indices = []
            for o in old:
                replace_indices.append(amp.indices.index(o))
            indices = deepcopy(amp.indices)
            for replace, n in zip(replace_indices, new):
                indices[replace] = n
            t_indices.append(indices)
        return t(t_indices)


class E(amplitude):
    def __init__(self, indices: list[str]) -> None:
        """Excitation operator class."""
        super().__init__(indices)

    def __str__(self) -> str:
        """String representation of E."""
        string = ""
        for a, i in zip(self.indices[::2], self.indices[1::2]):
            string += f"E_{{{a}{i}}} "
        return string

    def __add__(self, __value):
        """Addition of excitation operators."""
        assert type(__value) is E, "Excitation operators can only be added to other excitation operators"
        return E(self.indices + __value.indices)

    def __iter__(self) -> Generator[str, str]:
        """Iterate over excitation operators."""
        for a, i in zip(self.indices[::2], self.indices[1::2]):
            yield a, i

    def update_indices(self, old_indices: list[str], new_indices: list[str], translate=False):
        if not translate:
            for old, new in zip(old_indices, new_indices):
                assert old[0] == new[0], f"Old and new index must both be of the same type - here they were {old[0]} and {new[0]}."
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
        return E(indices)


class BRA(amplitude):
    def __init__(self, indices: list[str], is_excitation_vector: bool) -> None:
        """Bra class."""
        super().__init__(indices)
        self.left_excitation_vector = is_excitation_vector
        self.is_HF = self.indices == []

    def __str__(self) -> str:
        """String representation of bra."""
        if self.indices == []:
            return "<HF|"
        string = "<{}^{"
        for a in self.indices[::2]:
            string += a
        string += "}_{"
        for i in self.indices[1::2]:
            string += i
        string += "}|"
        return string

    def update_indices(self, old_indices: list[str], new_indices: list[str], translate=False):
        if not translate:
            for old, new in zip(old_indices, new_indices):
                assert old[0] == new[0], f"Old and new index must both be of the same type - here they were {old[0]} and {new[0]}."
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
        return BRA(indices, self.left_excitation_vector)


