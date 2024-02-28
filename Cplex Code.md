Cplex Code:

```
int T=...;  //number of Tasks
int R=...;  //number of Robots
int P=...;  //number of Positions

range t_range=1..T;
range p_range=1..T;
range r_range=1..R;

int d[t_range][t_range]=...;  //distance:Task t to t1
int D[r_range][t_range]=...;  //distance:Robot r to Task t
int h[t_range]=...;  //handling time of the rack associated with task t
int M[r_range]=...;  //endurance mileage of robot r

dvar boolean x[t_range][p_range][r_range];  //decision variable
dvar float S[p_range][r_range];
dvar float Cmax;

minimize Cmax;

subject to
{
forall (t in t_range)
  sum(r in r_range,p in p_range)
    x[t][p][r]==1;  //Constraint(1)
forall (p in p_range,r in r_range)
  sum(t in t_range)
    x[t][p][r]<=1;  //Constraint(2)
sum(t in t_range,r in r_range,p in p_range)
  x[t][p][r]==T;  //Constraint(3)
forall (p in (1..T-1),r in r_range)
  sum(t in t_range) x[t][p][r]>=sum(t in t_range) x[t][p+1][r];  //Constraint(4)
forall (r in r_range)
 S[1][r]>= sum(t in t_range) D[r][t]*x[t][1][r];  //Constraint(5)
forall (p in (2..T),r in r_range)
  S[p][r]>=S[p-1][r]+(sum(t in t_range) h[t]*x[t][p-1][r])+
  (sum(t in t_range,t1 in t_range) d[t][t1]*x[t][p][r]*x[t1][p-1][r]);  //Constraint(6)
forall (p in p_range,r in r_range)
  Cmax>=S[T][r]+sum(t in t_range) h[t]*x[t][T][r];  //Constraint(7)
forall (r in r_range)
  (sum(p in p_range,t in t_range) h[t]*x[t][p][r])+
  (sum(t in t_range) D[r][t]*x[t][1][r])+
  (sum(p in(1..T-1),t in t_range,t1 in t_range) (d[t][t1]*x[t][p][r]*x[t1][p+1][r]))<=M[r];  //Constraint(8)
forall (p in p_range,r in r_range)
  S[p][r]>=0; //Constraint(10)
Cmax>=0;  //Constraint(11)
}
```

