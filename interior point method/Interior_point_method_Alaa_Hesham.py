#!/usr/bin/env python
# coding: utf-8

# In[100]:


from sympy import *
import numpy as np


# #### Enter problem in the original form where A,c matrix only contain coefficient of design variables and s has same dmension as c.
# 

# In[101]:


def calculate_y(A_,c_,s_):

    A_=Matrix(A_).T
    c_=Matrix(c_)
    s_=Matrix(s_)
    y=A_.inv()*(c_-s_)    
    return y

A_=[[2, 1],[1 ,3]]
c_=[-30.0,-20.0]
s_=[1,1]
y=calculate_y(A_,c_,s_)
y



len(y)

#### Enter problem in standard form where A ,c matrix contain coeffiecient of slack variables where n is number of design variables   without slack variables , m number of constriants 


n = 2 # Number of design variables 
m = 2  # Number of constraints 
c= [-30.0,-20.0,0.0,0.0] # Coefficient of f 
A=[[2, 1, 1 , 0],[1 ,3, 0, 1]]
b=[8,8]
x=[2.5366,1.7561,1.1707,0.191]
s=[1,1,1,1]

sigma = 0.1
alpha = 0.3


# In[105]:


def attain_xys_length(x,y,s):
    x_l=len(x)
    y_l=len(y)
    s_l=len(s)
    return x_l,y_l,s_l


# In[106]:


x_l,y_l,s_l=attain_xys_length(x,y,s)


# In[ ]:





# In[107]:


def diag(x):
    l_x=len(x)
    X=np.zeros([l_x,l_x])
    for i in range (l_x):
        X[(i,i)]=x[i]
    return X 


# In[108]:


def introduce_inputs(A,b,c,x,X,s,S,y):
    A=Matrix(A)
    b=Matrix(b)
    c=Matrix(c)
    A_T=Matrix(A).T
    x=Matrix(x)
    X=Matrix(X)
    s=Matrix(s)
    S=Matrix(S)
    y=Matrix(y)
    return A,b,c,A_T,x,X,s,S,y


# In[109]:


def compose_Jacobian(m,n,A_T,A,S,X):
    
    # These twwo lines for making elements of A of float datatype
    A=np.array(A.tolist()).astype(np.float64) 
    A=Matrix(A)
    
    A_T=np.array(A_T.tolist()).astype(np.float64)
    A_T=Matrix(A_T)
    
    m_1_1=Matrix(np.zeros([m+n,m+n]))
    
    m_1_2=A_T
    m_1_3=Matrix(np.identity(m+n))
    
    m_2_1=A
    m_2_2=Matrix(np.zeros([m,m]))
    m_2_3=Matrix(np.zeros([m,m+n]))
    
    m_3_1=S
    m_3_2=Matrix(np.zeros([m+n,m]))
    m_3_3=X
    
    J=Matrix([[m_1_1,m_1_2,m_1_3],[m_2_1,m_2_2,m_2_3],[m_3_1,m_3_2,m_3_3]])
    
    return J

    
    


# In[110]:


def calculate_mu(x,s,n):    
    mu = (x.T*s)[0]/n
    return mu

def calculate_e(m,n):    
    
    e=Matrix(np.ones([m+n,1]))
    return e 


# In[111]:


def compose_F(A_T,y,s,c,A,x,b,X,S,e,sigma,mu):
    
    rc=[A_T*y+s-c]
   
    rb=[A*x-b]
    last_term=[X*S*e-(sigma+mu)*e]
    
    F = Matrix([rc,rb,last_term])
        
    return F


# In[112]:


def compose_xys(x,y,s):
    xys=Matrix([[x],[y],[s]])
    return xys
    


# In[113]:


def make_matrix_one_unit(R):
    R = BlockMatrix(R)
    R = R.as_explicit()
    return R


# In[114]:


def solver(J,F):
    
    delta_xys=J.inv()*(-1*F)
    
    return delta_xys
    


# In[115]:


def decompose_xys(xys,x_l,y_l,s_l):  # input to this function is sympy matrix, output is list of x, y, s
    
    x=list(xys[0:x_l])
    
    y=list(xys[x_l:x_l+y_l])
    
    s=list(xys[x_l+y_l:x_l+y_l+s_l])
    
        
    return x,y,s


# In[116]:


A_=[[2, 1],[1 ,3]]
c_=[-30.0,-20.0]
s_=[1,1]
y=calculate_y(A_,c_,s_)


# In[117]:



n = 2 # Number of design variables 
m = 2  # Number of constraints 
c= [-30.0,-20.0,0.0,0.0] # Coefficient of f 
A=[[2, 1, 1 , 0],[1 ,3, 0, 1]]
b=[8,8]
x=[2.5366,1.7561,1.1707,0.191]
s=[1,1,1,1]
sigma = 0.1
alpha = 0.3
print(x)


# In[118]:


X=diag(x)
S=diag(s)
A,b,c,A_T,x,X,s,S,y=introduce_inputs(A,b,c,x,X,s,S,y)
mu=calculate_mu(x,s,n)
e=calculate_e(m,n)

J=compose_Jacobian(m,n,A_T,A,S,X)

F=compose_F(A_T,y,s,c,A,x,b,X,S,e,sigma,mu)
xys=compose_xys(x,y,s)
print (x)
x_l,y_l,s_l=attain_xys_length(x,y,s)


# In[ ]:





# In[119]:


def unfold_update_xys(xys,alpha,delta_xys):
    
    xyy = []
   
    for t in xys:
        xyy.append(t)
        
    xys = Matrix(xyy)
    xys = xys+ alpha*delta_xys
    
    return  xys


# In[120]:


first_design_variable=[]
second_design_variable=[]
x_1=x[0]
x_2=x[1]
first_design_variable.append(x_1)
second_design_variable.append(x_2)

i=0
while (i<1):
    
    mu=calculate_mu(x,s,n)
    # make jacoobian one unit i.e. one matrix 
    J=BlockMatrix(J)
   
    F=BlockMatrix(F)    
    J = J.as_explicit()
    print (type(J))
    F = F.as_explicit()
    
    # solve for delta x,y,s
    delta_xys=J.inv()*(-1*F)
    
    # make xys one vector and update values of x,y,s 
    
    xys=unfold_update_xys(xys,alpha,delta_xys)
    
    # separate components of xys to allow them update jacobbian
    x,y,s=decompose_xys(xys,x_l,y_l,s_l)
    
    # prepare arguments for updated Jacobiaan 
    
    X=diag(x)
    S=diag(s)
    
    x_1=x[0]
    x_2=x[1]
    first_design_variable.append(x_1)
    second_design_variable.append(x_2)
    
    x=Matrix(x)
    s=Matrix(s)
    
    y=Matrix(y)
    J=compose_Jacobian(m,n,A_T,A,S,X)
    F=compose_F(A_T,y,s,c,A,x,b,X,S,e,sigma,mu)
    i+=1


# In[121]:


first_design_variable


# In[122]:


second_design_variable


# In[126]:


f_1=30*2.53660000000000+20*1.75610000000000
f_1


# In[127]:


f_2=30*4.72791780074959+20* 0.804674913878612


# In[128]:


f_2


# In[129]:


#we can see that using interior point method, we managed to maximize function to 157 


# In[ ]:




