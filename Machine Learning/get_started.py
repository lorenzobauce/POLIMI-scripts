# Assigning grades to students
grades = [20, 21, 25, 27, 18, 30, 17, 26, 28]
for i in grades:
    if i < 18:
        print("Fail")
    elif i >= 18 and i <= 21:
        print("Sufficient")
    elif i >= 22 and i <= 24:
        print("Fair")
    elif i >= 25 and i <= 27:
        print("Good")
    elif i >= 28 and i <= 30:
        print("Excellent")
