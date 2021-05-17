import qplot
import matplotlib.pyplot as plt

def plots_radiation_trapping():
    time = [5, 10, 25, 50]
    df_ne50 = "../results_radiation_trapping_ne_50_a0_2500/"
    df_ne100 = "../results_radiation_trapping_ne_100_a0_2500/"
    df_ne50_a0_150 = "../results_radiation_trapping_ne_50_a0_150/"
    for i in time:
        #qplot.density(df=df_ne50, clf=True, t=i)
        #plt.savefig("density_a0_2500_ne50_" + "t=" + str(i) + ".png")
            
        qplot.density(df=df_ne100, clf=True, t=i)
        plt.savefig("density_a0_150_ne50_" + "t=" + str(i) + ".png")

   # qplot.energy(df=df_ne50, clf=True)
   # plt.savefig("energy_a0_2500_ne50.png")
    qplot.energy(df=df_ne100, clf=True)
    plt.savefig("energy_a0_150_ne_50.png")

def main():
    plots_radiation_trapping()

if __name__ == "__main__":
    main()

