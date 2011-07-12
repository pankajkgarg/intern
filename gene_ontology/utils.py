"""Utility classes and routines that do not fit anywhere else"""

class bidict(object):
    """Bidirectional many-to-many mapping"""

    def __init__(self, items=None):
        object.__init__(self)

        self.left = {}
        self.right = {}

        if items is None:
            return

        if isinstance(items, bidict):
            # Copying an existing bidict
            self.left = dict(items.left)
            self.right = dict(items.right)
        elif isinstance(items, dict):
            # items assigns an element from the left dict to a
            # set of elements from the right dict
            for key, values in items.iteritems():
                self.add_left_multi(key, values)
        else:
            raise TypeError("items must be dict or bidict")

    def add_left(self, v1, v2):
        """Adds a pair of items `v1` and `v2` to the mapping s.t. `v1` as a
        left item is mapped to `v2` as a right item."""
        try:
            self.left[v1].add(v2)
        except KeyError:
            self.left[v1] = set([v2])
        try:
            self.right[v2].add(v1)
        except KeyError:
            self.right[v2] = set([v1])

    def add_right(self, v1, v2):
        """Adds a pair of items `v1` and `v2` to the mapping s.t. `v1` as a
        right item is mapped to `v2` as a left item."""
        return self.add_left(v2, v1)

    def add_left_multi(self, v1, v2s):
        """Associates multiple items in `v2s` to `v1` when `v1` is interpreted
        as a left item"""
        try:
            self.left[v1].update(v2s)
        except KeyError:
            self.left[v1] = set(v2s)
        for v2 in v2s: 
            try:
                self.right[v2].add(v1)
            except KeyError:
                self.right[v2] = set([v1])

    def add_right_multi(self, v2, v1s):
        """Associates multiple items in `v1s` to `v2` when `v2` is interpreted
        as a right item"""
        self.right[v2].update(v1s)
        for v1 in v1s:
            self.left[v1].add(v2)

    def get_left(self, v1, default = None):
        """Returns the items associated to `v1` when `v1` is looked up from the
        left dictionary"""
        return self.left.get(v1, default)

    def get_right(self, v1, default=None):
        """Returns the items associated to `v1` when `v1` is looked up from the
        right dictionary"""
        return self.right.get(v1, default)

    def len_left(self):
        """Returns the number of unique left items"""
        return len(self.left)

    def len_right(self):
        """Returns the number of unique right items"""
        return len(self.right)

    def iteritems_left(self):
        """Iterates over the left dictionary"""
        return self.left.iteritems()

    def iteritems_right(self):
        """Iterates over the right dictionary"""
        return self.right.iteritems()
