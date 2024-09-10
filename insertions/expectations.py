class Expectation:
	def __init__(self):
		self.expected = {}
		self.total = 0
		self.name = ""

	def get(self, minimum):
		"""Return the expected number and count at this minimum homology"""
		return self.expected[minimum], self.total


class ActinomycesExpMC(Expectation):
	"""Based on ANaesl, AVisc and AIsrael together"""
	def __init__(self):
		self.expected = {
				3: 120014,
				4: 34692,
				5: 8649,
				6: 2042,
				7: 496,
				8: 101,
				}

		# Based on balanced nucleotide sampling
# 		self.expected = {
# 				3: 132070,
# 				4: 38706,
# 				5: 9681,
# 				6: 2377,
# 				7: 552,
# 				8: 123,
# 				}

		self.total = 3000000
		# This was a simulation of random insertion positions, so we call it MC
		self.name = "Actinomyces MC"


class ANaeslExpMC(Expectation):
	"""Just based on ANaesl"""
	def __init__(self):
		self.expected = {
				3: 121210,
				4: 35356,
				5: 9103,
				6: 2154,
				7: 509,
				8: 116,
				}
		self.total = 3000000
		self.name = "A.Naeslundii MC"


class CodExp(Expectation):
	"""What we actually find in cod with -fas"""
	def __init__(self):
		self.expected = {
				3: 1312,
				4: 404,
				5: 123,
				6: 41,
				7: 8,
				}
		self.total = 144044
		# This was actual homology counts in Cod, no MC involved
		self.name = "Cod"


class HumanExp(Expectation):
	"""What we actually find in human with -fas"""
	def __init__(self):
		self.expected = {
				3: 7174,
				4: 2380,
				5: 629,
				6: 184,
				7: 76,
				}
		self.total = 235614
		self.name = "Human"


class PangolinExp(Expectation):
	def __init__(self):
		self.expected = {
				3: 3605,
				4: 1198,
				5: 433,
				6: 120,
				7: 0,	# We just didn't compute this
				}
		self.total = 75535
		self.name = "Pangolin"


class CodExpMC(Expectation):
	def __init__(self):
		self.expected = {
				3: 158642,
				4: 49952,
				5: 16163,
				6: 4644,
				7: 1345,
				8: 409,
				}
		self.total = 3000000
		self.name = "Cod MC"


class HumanExpMC(Expectation):
	def __init__(self):
		self.expected = {
				3: 163033,
				4: 50930,
				5: 15695,
				6: 4718,
				7: 1388,
				8: 387,
				}
		self.total = 3000000
		self.name = "Human MC"


class PangolinExpMC(Expectation):
	def __init__(self):
		self.expected = {
				3: 164421,
				4: 51529,
				5: 15999,
				6: 4688,
				7: 1406,
				8: 443,
				}
		self.total = 3000000
		self.name = "Pangolin MC"


class RabbitExpMC(Expectation):
	def __init__(self):
		self.expected = {
				3: 161346,
				4: 49873,
				5: 15187,
				6: 4412,
				7: 1212,
				8: 366,
				}
		self.total = 3000000
		self.name = "Rabbit MC"


class BatExpMC(Expectation):
	def __init__(self):
		self.expected = {
				3: 164513,
				4: 51508,
				5: 15769,
				6: 4587,
				7: 1338,
				8: 444,
				}
		self.total = 3000000
		self.name = "Bat MC"


class LizardExpMC(Expectation):
	def __init__(self):
		self.expected = {
				3: 161293,
				4: 50725,
				5: 15478,
				6: 4287,
				7: 1316,
				8: 429,
				}

		self.total = 3000000
		self.name = "Lizard MC"


class CodExpShuffle(Expectation):
	def __init__(self):
		self.expected = {
				3: 873,
				4: 236,
				5: 59,
				6: 17,
				}
		self.total = 32282
		self.name = "Cod Shuffle"


class HumanExpShuffle(Expectation):
	def __init__(self):
		self.expected = {
				3: 4069,
				4: 1190,
				5: 355,
				6: 97,
				}
		self.total = 85441
		self.name = "Human Shuffle"


class ANaeslExpShuffle(Expectation):
	def __init__(self):
		# Based on 500 iterations, reshuffling on each
		self.expected = {
				3: 1294,
				4: 331,
				5: 89,
				6: 20,
				}

		self.total = 42557
		self.name = "ANaesl Shuffle"


class AIsraelExpShuffle(Expectation):
	def __init__(self):
		self.expected = {
				3: 1838,
				4: 477,
				5: 132,
				6: 44,
				}
		self.total = 57083
		self.name = "AIsrael Shuffle"


class AViscExpShuffle(Expectation):
	def __init__(self):
		self.expected = {
				3: 1537,
				4: 386,
				5: 109,
				6: 20,
				}

		self.total = 49650
		self.name = "AVisc Shuffle"


class TreponemaExpShuffle(Expectation):
	def __init__(self):
		self.expected = {
				3: 1198,
				4: 402,
				5: 130,
				6: 35,
				}

		self.total = 18553
		self.name = "Treponema Shuffle"


class AActinomExpShuffle(Expectation):
	def __init__(self):
		self.expected = {
				3: 767,
				4: 236,
				5: 82,
				6: 27,
				}

		self.total = 13800
		self.name = "AActinom Shuffle"


class PorphyromonasExpShuffle(Expectation):
	def __init__(self):
		self.expected = {
				3: 1020,
				4: 326,
				5: 104,
				6: 27,
				}

		self.total = 20889
		self.name = "Porphyromonas Shuffle"


class TForsythExpShuffle(Expectation):
	def __init__(self):
		self.expected = {
				3: 1408,
				4: 388,
				5: 121,
				6: 35,
				}

		self.total = 26452
		self.name = "TForsyth Shuffle"
