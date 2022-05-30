# x = 3 + 4j
# print(abs(x))

# def flatten(mat):
#     out = []
#     for k in range(len(mat)):
#         for l in mat[k]:
#             out.append(l)
#     return out

# x = [[2], [3,4], [1], [5, 6, 7]]
# x = flatten(x)
# print(x)

x = {'100': 1, '110':-3, '111':-4, '000':10}
y = dict(sorted(x.items(), key=lambda item: item[0]))
print(y)