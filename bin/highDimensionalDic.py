from collections import defaultdict

def nestedDic():
    return defaultdict(nestedDic)

def OneDepDic():
    return defaultdict(list)

def TwoDepDic():
    return defaultdict(OneDepDic)

def ThreeDepDic():
    return defaultdict(TwoDepDic)

def FourDepDic():
    return defaultdict(ThreeDepDic)

def FiveDepDic():
    return defaultdict(FourDepDic)