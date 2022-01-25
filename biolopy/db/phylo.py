import re

bep = "((osat,lper),(bdis,(atau,hvul)))"
pacmad = "((sbic,zmay),(sita,phal))"
poaceae = f"({bep},{pacmad})"
monocot = f"(({poaceae},macu),drot)"


def extract_species(tree: str):
    return re.findall(r"[^(), ]+", tree)
