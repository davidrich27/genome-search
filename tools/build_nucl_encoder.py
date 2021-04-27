
nucl_dict = {
   "A": 0,
   "C": 1,
   "G": 2,
   "T": 3,

   "a": 0,
   "c": 1,
   "g": 2,
   "t": 3
}

print("int nucl_encoder = { ", end="")

for i in range(256):
   c = chr(i)
   if c in nucl_dict.keys():
      print("{}".format(nucl_dict[c]), end="")
   else:
      print("-1", end="")
   if i < 255:
      print(", ", end="")

print(" };")