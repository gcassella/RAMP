# MCRAMP
Monte Carlo Raytracing Achieved via Massive Parallelisation

![example](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAqMAAAHcCAYAAADiA6PhAAAABHNCSVQICAgIfAhkiAAAAAlwSFlzAAAPYQAAD2EBqD+naQAAIABJREFUeJzs3Xt8U/X9x/F3gFLuLddCuVYQRBRQJohX1E5AJ+JlU9wFt6kTwXmdiAooc2PTOZ2OyXQqu3h3ispPUUFAEQQVEZGLoIjcytW2QGkL7ff3x3enOUmTNmmTnqR5PR+PPJImJ+k3p2nyzud7OT5jjBEAAADggQZeNwAAAACpizAKAAAAzxBGAQAA4BnCKAAAADxDGAUAAIBnCKMAAADwDGEUAAAAniGMAgAAwDOEUQAAAHiGMAoAAADPEEYBAADgGcIoAAAAPEMYBQAAgGcIowAAAPAMYRQAAACeIYwCAADAM4RRAAAAeIYwCgAAAM8QRgEAAOAZwigAAAA8QxgFAACAZwijAAAA8AxhFAAAAJ4hjAIAAMAzhFEAAAB4hjAKAAAAzxBGAQAA4BnCKAAAADxDGAUAAIBnCKMAAADwDGEUAAAAniGMAgAAwDOEUQAAAHiGMAoAAADPEEYBAADgGcIoAAAAPEMYBQAAgGcIowAAAPAMYRQAAACeIYwCAADAM4RRAAAAeIYwCgAAAM8QRgEAAOAZwigAAAA8QxgFAACAZwijAAAA8AxhFAAAAJ4hjAIAAMAzhFEAAAB4hjAKAAAAzxBGAQAA4BnCKAAAADxDGAUAAIBnCKMAAADwDGEUAAAAniGMAgAAwDOEUQAAAHiGMAoAAADPEEYBAADgGcIoAAAAPEMYBQAAgGcIowAAAPAMYRQAAACeIYwCAADAM428bgAA1EfGSMXF0nffSfn50v790oED9nTwoD0vKpJKSgJPpaVSWZlUXu4/lZVJPp/UsKHUoIE9b9hQatRISk+XmjTxnzdpIjVvLrVsKbVo4T/PzJRat5aaNrWPBQCJgjAKABEqLZV27JC2b5e2bZPy8qRdu+xp9257vmePDaDffWe3D8fnk5o1syHSfWrc2B86neDphMeyMn9QLSuTjhyxAba42H9eXGxvCyc93YbSNm2ktm2lrKzAU8eOUpcuUteuUrt2BFcA8eczxhivGwEAieDQIWnjRmnTJmnzZumbb/ynrVtt2HRr3Fjq0EFq395/3q6dDXqtW/urkZmZUqtWtkLZooWtXMarQmmMDcFOJXb/fnsqKJD27bMh2TnfvVvaudN/2rPHBl3383OC6VFHBZ569bLPFQBqizAKIKUYY6uaX3xhT19+aU8bNtjA6UhPl7p3l3r0sOddu0qdO9tTdrY9b926flUOy8ps4N62ze6LLVvs+ebNNqB//bUNrI62baW+faVjjrGnfv2kAQNsdbU+7RcA8UUYBVBvFRVJn38uffqpPX3+ubRmja0SSrY62bu3dPTR9uRc7tnTVjobMMWzksJCG0y//FJaty7wVFRkt2nf3obSAQOkk06SBg+2oZ6ACiAUwiiAeuHwYRs2ly2zp48+sgGpvNyOuzz2WKl/f1u9c049etjbUHvl5TakrlolffaZPa1caYc4SLaKOniwdMop0umn28tNm3raZAAJgjAKICnt3y8tWSK99570/vvSxx/bMZ+NGtmK3ODB0gkn2NNxx9lZ5qh7u3fbLwbLl9svCUuX2sp048a2ajpsmHTuudLJJ9vrAKQewiiApHL55TZ87txpxzh26GArbaecIg0ZIp14YmJX3MrK/Es7HTwY2dJO5eV2rKuzpJP75Czn5D45Szs1SsD1UsrKbAX7/fftF4kFC6S9e+3ErrPOkkaMkEaNshOnAKQGwiiApDJkiPTVV9LvfiedeabUp483YxGPHPEv6xR82rfPri3qLPGUn2/HWh44YKu3daVJE/9ao5mZdpa/+9SunX85J+e8TZu6HStbXm7H8779tj0tXmz37UknSaNHSxdfbCdHAai/CKMAksoNN0jz50urV8fvdxhjq3WbNvmXdtq82c4yd055eYHLIEk2+LVv71/ayb28U6tW/gXomze3582a+Resd68z2qiRf51R5yT51xl1KqfudUadtUYPHfIvru8s67R/vw3E+/b5T3v32pnxBw8GPoe0NLtSQNeu/lP37nY5p5497eW0tPjt+/x86Y03pNmzpTfftM/jxBOlH//YVsWzs+P3uwF4gzAKIKlMmiQ9/7xdZqi2vvtOWrvWLuvkLPH05Ze28uoOaS1b2hDWpYt/eSdniaesLP8ao82a1b5Nde3AgcC1Rrdvt0s6uU9bt9rgK9mhAd262VUH+va1p2OPteexXne0uNgG06eflubMsZPUzjtPuvZaaeRIJp8B9QVhFEBS+e1vpRkzbGUyUiUldk3Rzz6zFVXntH27f5suXezSTr172wXdc3LsbPsePerfeqLROnLE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1pjOnY3p2tWY5cvDP+bw4fY+l15qzzMy/G384gv/3+jPfzbm5z+370F33WXfFyVjcnLs45SWGtOokf++XbsGPt8HHvD/zrPOstf94x9VP9/evf2v2aFD/Y914on29v/+139d69b+7as65eaGv23AgMrXPfhgYJvmzvXvJ/djtWnj3+a99/zXf/CBvc557P/7v/DP9/TT7TbTpgVeX17uv1xYaMxRR9ntbrjBmLKy8I93662hn+err/q3KSryv94vu8xet2+fMccdF/q+zz1nt7nlltC3d+pkzJNPGnPokP93bN9uzAknVN62TRtjpk+3n8tun39uzF//aj+jx4yx23boEPo5Ll8euh27doXe3nm9O//7CxYYs2VL4D42xpj160O3WTKme3djjhyp/Nh33x243ahR4dvheOstY3r18t/nyisrtyUeapvFfMZUl+Xjp7S0VM2aNdNLL72k0aNHV1w/duxY5efn69VXX610n5KSEpW4vpoXFhaqa9euKigoUCvn622cLF1aebDyzp12qYzmzaXHH7ffNjt0sGX78nL7bX3OnMD7dO8unXRS6LcvqXbX5+VVrg788If+y+6/tjGV1zP70Y9qtm+QOjZssDOlHRkZgdUct+++k+bNC7zukkv8a0TG2pw5gd1SHTr4q0qbN9tqzYgR0tlnS7fdZitizrI4VXn5ZVsV3LTJduW9/76tkg4YYKtVHTva9wGnG/38821FZ8MGacoUW0F6+WU7ZvLQIVvh69vXvmecdpqtFo0a5R+T6ubs77Zt7cQKZ9xfWpqtvvbta6tOBQW269WplFXHGOnFF+3lvXvtZAl3dWbMGDtO88orbbtatPD38jRsKJ11lt1/TvUm+DEXLfJXBx1FRbZb0pnZHcoJJ/iXBjr9dNvNL9nK78kn+7vDH3vM3+18223S/ffbv4NTmQtWUiI5HynffmsrXs76yrfeau+/f7+tIh45Ypctu+wy6fvf9z/G00/biVT//Kf9uUMHW6Hs398/nvPkk21V91//8t/vtdfsc3JmR59wgv0caNjQ/x581VW2a9aZkHXVVfb14ezXY47xVwnPO89Wa0tK7GvLGV4Q6jm//779PHLGSYfy6quS8xHco4etwKel2f2Rn297Nw4c8Fdju3f396o1a2Yrw82b+x9v2TK7HyS7RJG7Cn711dLzz9sehKws2yvQrJn9vFy82O7HPn1sBfO22/yV9qZNbdd4q1Z2BYnt2+3wh9des/vgnnv8vXtpafYztm1b285Q42EHDw5/BKW//S3wf+HEE8OvEDBnjh2CEaxjR1uxbdHC/n+//77dh82a2X3dvr0dGrFjh62oO5PJgv3rX3a1DHePQ58+9vHbtrWP7WSNoiL7ujDG3v7EE+EfN9YKCwuVkZFR8ywW02gcpW3bthlJZsmSJQHX/+Y3vzGDBw8OeZ+pU6caSZVOdVEZff758N+AJ02y34YaNgx9e+fO4b9VxuN0003G+Hx19/s4carqdOmltsJRF7+rRw9jJk4Mf/uvfmWrfM2aRfe4nTvbCsOf/lT5tiFD7HvEL35R+bb/+11BdGsAABl/SURBVD9jVq40pmXLwOtHjrRVqIsuiuz3//KXxhQXG9O3b+D1d9xhzD//aUyDBjXbX1lZ9nnt3Gmfo3P9H/5gf1///oHbh6vuBO+rgweNmT3bf51TMT5wwF8l8/mMueoqW6VytnNXPJcuNeb44435y1/sz4895t/OXZHfssWYdu0ie76nnmqf75Ejxlx/vTEnn2zMN9/4H2vOHFupzcuzP0+fbqvSZ55pq7UlJcacf759rN/+1v/7Tz3VVljXrrWVuZ/+1P5NfvYz/++76y5j0tJCt8upNm/caN+/3W0yxlYDgyu1LVva/Vxb5eW2CpeZWf3+u/pqe5/du23FNVxl9oEHAj+DMjONWbXK3lZUZMzixcbs2VN1u/Lzjfn1r0O3KzPT7mtHaan9nTk5lbdt3Nj2Ptx6q30vGjjQmFmzqv7dc+bY/Z2WZsxLL1W97Z49xrzwgr93IdzptNOM2bEj8L6RVC2Li6t/bPdpzBj796lLSV0Z3b59uzp37qwlS5ZoqDOKWtJtt92mRYsWaVmIry1eVkY3bAg9yLt3b1sZ8vnsETC++MIODG/Y0J43bmwrDO3a2W9RW7b4l+Hw+UKfanNb27b227zTlmDuJUDS0+23tDVr4rNECeqnzEy7MPOCBVWPd5Ts6/+SS+w399deC1zvM9Z8PjvOrFcvOxHkq6/s27MjPd22pW1bW+V8663KSzCFc/bZ/grUp5/aSRbG2P/xkSNt1dMYWwFx/peysmzPhM9nKznz5tkKZosWdnHrjAy7X+bOtVWncPulWTPb69K0qa0ovf22rTqnpdnHb97cVs0++CD6Qw6feaat7En2sV95xbZr1Cj73IqLbbu3brVjLTt3tlWizz+3VV73WsWSbdOoUf4JnWvW2KrgySf7Z74XFtpq3LHH2jGC+/fb96sDB2zvTFWV85077b4aPDjwvWzbNls5rurv2by5fXznaFmRKi21bXLaZYytGGZmBrYhWEmJfyKYY8cOWzUuKLD7rmFDW2F3KonVtePNN+1jNGlix7g6k8JioajI/q337rX7sUUL+xybNLEV8MxMO8mqqufstmmT/V/Iz7cV/VDjgCNRXm4nPDmvOcn+nwdPKpLs32bjRvsaPXjQVlKHDfOvlBCNw4dtlbyqqnKwwkL7+ty+3VZ3jxyxp7ZtbY+JuychWiUldp/m5dnT3r32b+FkjYYN7fjd3NzI/0axUtvKaNJ10werdWkYAAAANVbbLObpbPrGjRtr0KBBmj9/fsV15eXlmj9/fkClFAAAAPVTLQrGsXHzzTdr7Nix+t73vqfBgwfroYce0sGDB/Xzn//c66YBAAAgzjwPo5dddpl2796tKVOmKC8vTwMHDtTcuXOV5SxGBwAAgHrL0zGjscCYUQAAAO8k9ZhRAAAApDbPu+lryynsFsbzGJsAAAAIyclgNe1sT/owuv9/B6Pt6j44LgAAAOrU/v37lZGREfX9kn7MaHl5ubZv366WLVvKV9ervCYR5+AAW7ZsYWxthNhn0WOf1Qz7LXrss+ixz6LHPouMMUb79+9Xdna2GjSIfgRo0ldGGzRooC5dunjdjKTRqlUr/qGixD6LHvusZthv0WOfRY99Fj32WfVqUhF1MIEJAAAAniGMAgAAwDMN77777ru9bgTqRsOGDTVs2DA1apT0ozPqDPsseuyzmmG/RY99Fj32WfTYZ/GX9BOYAAAAkLzopgcAAIBnCKMAAADwDGEUAAAAniGMAgAAwDOE0Xrum2++0S9/+Uvl5OSoadOm6tmzp6ZOnarS0tKA7Xw+X6XTc88951GrE8OMGTPUo0cPNWnSREOGDNHy5cu9blLCmD59uk466SS1bNlSHTp00OjRo7V+/fqAba688spKr6kRI0Z41GLv3X333ZX2xzHHHFNxuzFGU6ZMUadOndS0aVPl5uZqw4YNHrbYez169Aj53jR+/HhJvMYk6b333tMFF1yg7Oxs+Xw+zZ49O+D2SF5XxcXFGj9+vNq2basWLVrokksu0c6dO+vyadS5qvbb4cOHNXHiRB1//PFq3ry5srOz9bOf/Uzbt28PeIxhw4ZVev1de+21df1U6gXCaD23bt06lZeX6+9//7u++OILPfjgg5o5c6buuOOOSts+9dRT2rFjR8Vp9OjRHrQ4MTz//PO6+eabNXXqVK1YsUIDBgzQ8OHDtWvXLq+blhAWLVqk8ePH68MPP9Q777yjw4cP69xzz9XBgwcDthsxYkTAa+rZZ5/1qMWJoV+/fgH7Y/HixRW33XfffXr44Yc1c+ZMLVu2TM2bN9fw4cNVXFzsYYu99dFHHwXsr3feeUeS9MMf/rBim1R/jR08eFADBgzQjBkzQt4eyevqpptu0uuvv64XX3xRixYt0vbt23XxxRfX1VPwRFX7raioSCtWrNDkyZO1YsUKvfzyy1q/fr1GjRpVadurr7464PV333331UXz6x+DlHPfffeZnJycgOskmVdeecWjFiWewYMHm/Hjx1f8XFZWZrKzs8306dM9bFXi2rVrl5FkFi1aVHHd2LFjzYUXXuhhqxLL1KlTzYABA0LeVl5ebjp27Gjuv//+iuvy8/NNenq6efbZZ+uqiQnvhhtuMD179jTl5eXGGF5jwYLfxyN5XeXn55u0tDTz4osvVmyzdu1aI8ksXbq07hrvoUg+/5YvX24kmc2bN1dcd+aZZ5obbrgh3s1LCVRGU1BBQYHatGlT6frx48erXbt2Gjx4sJ588kmZFF2CtrS0VJ988olyc3MrrmvQoIFyc3O1dOlSD1uWuAoKCiSp0utq4cKF6tChg/r06aNx48Zp7969XjQvYWzYsEHZ2dk66qij9OMf/1jffvutJGnTpk3Ky8sLeM1lZGRoyJAhvOb+p7S0VP/5z3/0i1/8Qj6fr+J6XmPhRfK6+uSTT3T48OGAbY455hh169aN155LQUGBfD6fMjMzA65/+umn1a5dOx133HGaNGmSioqKPGphcuNwAilm48aNeuSRR/SnP/0p4Ppp06bp7LPPVrNmzfT222/ruuuu04EDB/TrX//ao5Z6Z8+ePSorK1NWVlbA9VlZWVq3bp1HrUpc5eXluvHGG3XqqafquOOOq7h+xIgRuvjii5WTk6OvvvpKd9xxh0aOHKmlS5eqYcOGHrbYG0OGDNGsWbPUp08f7dixQ/fcc49OP/10rV69Wnl5eZIU8jXn3JbqZs+erfz8fF155ZUV1/Eaq1okr6u8vDw1bty4UsjitedXXFysiRMnasyYMWrVqlXF9VdccYW6d++u7OxsrVq1ShMnTtT69ev18ssve9ja5EQYTVK33367/vjHP1a5zdq1awMmSGzbtk0jRozQD3/4Q1199dUB206ePLni8gknnKCioiLdf//9KRlGEZ3x48dr9erVAeMfJenyyy+vuHz88cerf//+6tmzpxYuXKhzzjmnrpvpuZEjR1Zc7t+/v4YMGaLu3bvrhRdeUN++fT1sWXJ44oknNHLkSGVnZ1dcx2sM8Xb48GH96Ec/kjFGjz76aMBt11xzTcXl448/XtnZ2Tr77LP11VdfqWfPnnXd1KRGN32SuuWWW7R27doqT0cddVTF9tu3b9dZZ52lU045RY899li1jz948GBt3bpVJSUl8XwaCaldu3Zq2LBhpdmkO3fuVMeOHT1qVWKaMGGC5syZowULFqhLly5VbnvUUUepXbt22rhxYx21LrFlZmaqd+/e2rhxY8XritdcaJs3b9a8efN01VVXVbkdr7FAkbyuOnbsqNLSUuXn54fdJlU5QXTz5s165513AqqioQwePFiSeP3VAGE0SbVv317HHHNMlafGjRtLshXRYcOGadCgQXrqqafUoEH1f/aVK1eqdevWSk9Pj/dTSTiNGzfWoEGDNH/+/IrrysvLNX/+fA0dOtTDliUOY4wmTJigV155Re+++65ycnKqvc/WrVu1d+9ederUqQ5amPgOHDigjRs3qlOnTsrJyVHHjh0DXnOFhYVatmwZrznZlT46dOig888/v8rteI0FiuR1NWjQIKWlpQVss379en377bcp/dpzguiGDRs0b948tW3bttr7rFy5UpJ4/dWExxOoEGdbt241vXr1Muecc47ZunWr2bFjR8XJ8dprr5nHH3/cfP7552bDhg3mb3/7m2nWrJmZMmWKhy331nPPPWfS09PNrFmzzJo1a8w111xjMjMzTV5entdNSwjjxo0zGRkZZuHChQGvqaKiImOMMfv37ze33nqrWbp0qdm0aZOZN2+eOfHEE83RRx9tiouLPW69N2655RazcOFCs2nTJvPBBx+Y3Nxc065dO7Nr1y5jjDF/+MMfTGZmpnn11VfNqlWrzIUXXmhycnLMoUOHPG65t8rKyky3bt3MxIkTA67nNWbt37/ffPrpp+bTTz81ksyf//xn8+mnn1bM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sbFRREY++CQy8gGmq+cgIpKZmSlms1m5/u233+TBBx8UtVotRqNRdu3aJTNmzJAdO3aM+p0dP35cYmNjRaPRiMlkks8//1yMRqPToScA8sUXXzi18/b2lvLycuW6urpaQkNDRa1Wy+LFi6WsrEwAyF9//SUiww8wiYgcOHBAEhISRKvVisFgkLi4OKd/B7jawYMHJTAwcNR1XDE4OCiZmZliMBjEx8dH1q1bN2z/xtuTiayHiCYHlcg1f/5GREQu09nZidmzZ+PgwYNYunTpvzr2li1bsGPHDly4cMEl/a1duxZ2ux1lZWUu6e96uXo9RHRr4GN6IiIX+vrrr9Hb24uoqChYrVasX78ewcHByu8nb6aysjIsWLAAfn5+sFgsKCoqwksvveSy/iMjIxEfH++y/sZzs9dDRLcGJqNERC40ODiIV199FWfPnoWXlxcSEhKwa9euYSfhb4b29nYUFBTgzz//RFBQELKzs5Gbm+uy/levXu2yvibiZq+HiG4NfExPRERERG7Dd9MTERERkdswGSUiIiIit2EySkRERERuw2SUiIiIiNyGySgRERERuQ2TUSIiIiJyGyajREREROQ2TEaJiIiIyG2YjBIRERGR2zAZJSIiIiK3YTJKRERERG7DZJSIiIiI3IbJKBERERG5DZNRIiIiInKb/wCfrZhMjTG1tAAAAABJRU5ErkJggg==)
*Above: visualization of instrument composed of guide, spherical powder sample, and detector. Below: Resulting histogram in the detector of the above instrument.*