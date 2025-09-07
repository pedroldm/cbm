import re

class Validator:
    lines: int
    columns: int
    binary_matrix: list[list[int]]

    def __init__(self, instance_path: str):
        with open(instance_path) as f:
            self.lines, self.columns = map(int, f.readline().split())
            self.binary_matrix = [[0 for _ in range(self.columns)] for _ in range(self.lines)]

            for i in range(self.lines):
                values = list(map(int, f.readline().split()))
                edges = values[1:]
                for e in edges:
                    self.binary_matrix[i][e - 1] = 1

    def validate(self, perm: list[int], cost: int):
        v1 = self.contains_all_columns(perm)
        v2 = self.count_one_blocks(perm)
        if not (v1 and (cost == v2)):
            raise Exception(f"Invalid solution: {v1} or {v2} != {cost}")
        return True

    def count_one_blocks(self, perm: list[int]) -> int:
        perm_matrix = [[row[j] for j in perm] for row in self.binary_matrix]
        bitstring = "-".join("".join(str(cell) for cell in row) for row in perm_matrix)
        blocks = re.findall(r'1+', bitstring)
        return len(blocks)

    def contains_all_columns(self, perm: list[int]) -> bool:
        return set(perm) == set(range(self.columns))