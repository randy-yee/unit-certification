read("ExperimentFunctions.py");
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
print("Experiment file 1-1");
{
sigstring = "1-1";

b_ranges =
[
    1,10,1;          \\#max 100, likely range [1,10]
    1,10,1;          \\#max 100, likely range [1,10,]
    1,10,1;          \\#3 max 100, likely range [1,10]
    30,60,3;          \\#max 400, likely range [30,70]
    5,35,3;           \\#max 200, likely range [5,35]
    4,25,2;           \\#6 max 200, likely range [1,30]
    15,60,3;           \\#max 400, likely range [60,160]
    60,160,8;          \\#max 800, likely range [60,160]
    40,120,10;         \\#9 max 400, likely range [40,120]
    300, 500, 15;       \\#max1500 , likely range [300,500]
    100, 300, 25;       \\#max600 , likely range [100,300]
    100, 400, 30;       \\#12 max600 , likely range [100,400]
    300, 700, 20;      \\#max 3300, likely range [600,800]
    300, 900, 50;     \\#max2600 likely range [500,900]
    350, 700, 50;     \\#15 ceil(i/3)== 5
    150,10000,20;      \\#ceil(i/3)== 6
    150,10000,20;      \\#ceil(i/3)== 6
    150,10000,20       \\#ceil(i/3)== 6
];
}
