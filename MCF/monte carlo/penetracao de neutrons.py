import numpy as np
import numpy.random as rd

ab = input("Item a ou b? ")

if ab == "a":


	pa = 1/3
	pe = 1 - pa
	t = 1

	lbd = [0.01, 0.05, 0.1, 1, 10]


	for lbd in lbd:

		A = 0
		R = 0
		T = 0


		for i in range(1000):

			i += 1

			z = 0
			
			while True:
				rand = rd.random()

				if rand < pa:
					A += 1
					break
				
				else:
			
					cos_theta = 1 - 2 * rd.random()
					l = -lbd * np.log(rd.random())
					z += l * cos_theta
					
					if z > t:
						T += 1
						break
					elif z < 0:
						R += 1
						break

		print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")      
		print("Lambda: ", lbd)
		print("Absorvidos: ", A)
		print("Refletidos: ", R)
		print("Transmitidos: ", T)	
		print("Probabilidade de absorção: ", A/1000)
		print("Probabilidade de reflexão: ", R/1000)
		print("Probabilidade de transmissão: ", T/1000)
		print("Probabilidade de absorção + reflexão + transmissão: ", (A+R+T)/1000)
		print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  
if ab == "b":
	
	pa = 0.5
	pe = 1 - pa
	t = 1
	lbd = 0.5

	A = 0
	R = 0
	T = 0


	for i in range(1000):

		i += 1

		z = 0
		
		while True:
			rand = rd.random()

			if rand < pa:
				A += 1
				break
			
			else:
		
				cos_theta = 1 - 2 * rd.random()
				l = -lbd * np.log(rd.random())
				z += l * cos_theta
				
				if z > t:
					T += 1
					break
				elif z < 0:
					R += 1
					break

	print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")      
	print("Lambda: ", lbd)
	print("Absorvidos: ", A)
	print("Refletidos: ", R)
	print("Transmitidos: ", T)	
	print("Probabilidade de absorção: ", A/1000)
	print("Probabilidade de reflexão: ", R/1000)
	print("Probabilidade de transmissão: ", T/1000)
	print("Probabilidade de absorção + reflexão + transmissão: ", (A+R+T)/1000)
	print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

	
else:
	print("Opção inválida. Por favor, escolha 'a' ou 'b'.")