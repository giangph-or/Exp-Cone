import numpy as np
nProducts = 200
nCustomers = [5,10,15,20]
nGroups = [5,10,20]
groupCapacity = [5,10,20]
id = [2,3]

for c in range (len(nCustomers)):
  for g in range (len(nGroups)):
    for gc in range (len(groupCapacity)):
      for idx in range (len(id)):
        productGroup = []
        revenue = np.random.uniform(low=1, high=3, size=(nCustomers[c],nProducts))
        utility = np.random.uniform(low=0, high=1, size=(nCustomers[c],nProducts))
        #fraction = np.random.uniform(low=0, high=1, size=(nCustomers))
        cost = np.random.uniform(low=1, high=2, size=(nProducts))

        for _ in range (0, nProducts):
          productGroup.append(np.random.randint(0, 2, size=(nGroups[g],)))

        f = open("C://Users//giang//source//repos//giangph-or//CPAOOverlap//x64//Release//AO_data//" + str(nProducts) + "_" + str(nCustomers[c]) + "_overlap_" + str(nGroups[g]) + "_" + str(groupCapacity[gc]) + "_" + str(id[idx]) + ".dat", "w")
        print(str(nProducts) + "_" + str(nCustomers[c]) + "_overlap_" + str(nGroups[g]) + "_" + str(groupCapacity[gc]) + "_" + str(id[idx]) + ".dat")
        # General capacity
        #f = open(str(nProducts) + "_" + str(nCustomers) + "_generalized_1_3.2" + ".dat", "w")
        f.write(str(nProducts))
        f.write("\n")
        f.write(str(nCustomers[c]))
        f.write("\n")
        f.write(str(nGroups[g]))
        f.write("\n")
        f.write(str(groupCapacity[gc]))
        f.write("\n")

        for i in range (0, nCustomers[c]):
          for j in range (0, nProducts):
            f.write(str(utility[i][j]) + " ")
          f.write("\n")

        for i in range (0, nCustomers[c]):
          for j in range (0, nProducts):
            f.write(str(revenue[i][j]) + " ")
          f.write("\n")

        # for j in range (0, nCustomers):
        #   f.write(str(fraction[j]) + " ")
        # f.write("\n")

        for j in range (0, nProducts):
          f.write(str(cost[j]) + " ")
        f.write("\n")

        for i in range (0, nProducts):
          for j in range (0, nGroups[g]):
            f.write(str(productGroup[i][j]) + " ")
          f.write("\n")

        # min = 10
        # max = 0

        # for j in range (0, nProducts):
        #   if revenue[j] < min:
        #     min = revenue[j]
        #   if revenue[j] > max:
        #     max = revenue[j]
        # f.write(str(min/max))

        f.close()