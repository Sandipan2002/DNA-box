'''
Author: Sandipan Chowdhury
0: 3' - 5'
1: 5' - 3'
'''
class NA():
	def __init__(self,seq,orient = 1):
		"seq=>sequence of neucleotides. orient=> 0 or 1(default)[0: (3'-5'), 1: (5'-3')"
		self.seq = seq
		self.orient = orient
		self._type = self.type_NA(self.seq)
		self.reference_create()
		self.translationDictionary = {"UUU": "phe", "UUC": "phe", "UUA": "leu", "UUG": "leu", "UCU": "ser", "UCC": "ser", "UCA": "ser", "UCG": "ser", "UAU": "tyr", "UAC": "tyr", "UAA": "end", "UAG": "end", "UGU": "cys", "UGC": "cys", "UGA": "end", "UGG": "trp", "CUU": "leu", "CUC": "leu", "CUA": "leu", "CUG": "leu", "CCU": "pro", "CCC": "pro", "CCA": "pro", "CCG": "pro", "CAU": "his", "CAC": "his", "CAA": "gln", "CAG": "gln", "CGU": "arg", "CGC": "arg", "CGA": "arg", "CGG": "arg", "AUU": "ile", "AUC": "ile", "AUA": "ile", "AUG": "met", "ACU": "thr", "ACC": "thr", "ACA": "thr", "ACG": "thr", "AAU": "asn", "AAC": "asn", "AAA": "lys", "AAG": "lys", "AGU": "ser", "AGC": "ser", "AGA": "arg", "AGG": "arg", "GUU": "val", "GUC": "val", "GUA": "val", "GUG": "val", "GCU": "ala", "GCC": "ala", "GCA": "ala", "GCG": "ala", "GAU": "asp", "GAC": "asp", "GAA": "glu", "GAG": "glu", "GGU": "gly", "GGC": "gly", "GGA": "gly", "GGG": "gly"}

	def getStrand(self,orient):
		"orient = 1 or 0, int. returns the strand of certain orient."
		if orient == 1:
			return self.seq
		if orient == 0:
			return self.complementStrandSyn(self.seq, res = self._type)
	def reference_create(self):
		"internal function to universalize"
		if self.orient == 1:
			pass
		elif self.orient == 0:
			self.seq = self.complementStrandSyn(self.seq,res = self._type)
			self.orient = 1
	def complementStrandSyn(self,strand, res = "DNA"):
		"res can be DNA(default) or RNA, returns Complementary Strand"
		if res == "DNA":
			mstrand = strand.replace("T","S")
			mstrand = mstrand.replace("C","B")
			mstrand = mstrand.replace("A","T")
			mstrand = mstrand.replace("G","C")
			mstrand = mstrand.replace("S","A")
			mstrand = mstrand.replace("B","G")
			return mstrand
		elif res == "RNA":
			mstrand = strand.replace("U","S")
			mstrand = mstrand.replace("C","B")
			mstrand = mstrand.replace("A","U")
			mstrand = mstrand.replace("G","C")
			mstrand = mstrand.replace("S","A")
			mstrand = mstrand.replace("B","G")
			return mstrand
	def type_NA(self,seq):
		"detects the type of NA. returns DNA or RNA."
		if "T" in seq:
			return "DNA"
		elif "U" in seq:
			return "RNA"
	def translate(self,seq):
		"takes a sequence as arg and returns translated hnRNA"
		x = seq
		if "AUG" in x:
			i = x.index("AUG")
			list_of_codons=[]
			while i+3<=len(x):
				list_of_codons+=[x[i:i+3]]
				i+=3
			peptide_seq = []
			for c in list_of_codons:
				peptide_seq+=[self.translationDictionary[c]]
			if "end" in peptide_seq:
				resonable_peptide = peptide_seq[0:peptide_seq.index("end")]
			else:
				resonable_peptide = peptide_seq
			return "-".join(resonable_peptide)
		else:
			return None
	def transcrypt(self,seq):
		return seq.replace("T","U")
	def get_self_translated(self,ori=1):
		"Takes ori. Returns the translated copy of the ori."
		return self.translate(self.get_self_transcrypted(ori))
	def get_self_transcrypted(self,ori=1):
		"takes ori and returns transcrypted copy."
		return self.transcrypt(self.getStrand(ori))


