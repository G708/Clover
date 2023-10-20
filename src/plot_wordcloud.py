import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from wordcloud import WordCloud
from sklearn import preprocessing


def asc4rank(score):
    """
    Return True if the score is in ascending order for ranking.
    This rank will reflect the word size in the word cloud.
    So, FDR and Glint will be False, the smaller the value, the larger the word size.
    and the others will be True, the larger the value, the larger the word size.

    Parameters
    ----------
    score : str
        A score name.

    Returns
    -------
    bool
        True if the score is in ascending order for ranking.

    """
    asc = ['FDR', 'Glint']
    if score in asc:
        return False
    else:
        return True

class FreqColorFunc(object):
    """
    A class that represents a color function for a word cloud based on the normalized score of each word.
    """
 
    def __init__(self, freq, colormap):
        """
            Parameters
           ----------
           freq : dict(str -> float)
           colormap : matplotlib.colors.Colormap 
 
        """
        self.m=255
        self.cm=colormap
        self.f=freq
        
    def c(self,x):
        """
        Returns an RGB color tuple for a given value.
        
        Parameters
        ----------
        x : float
            A value between 0 and 1.

        Returns
        -------
        tuple
            An RGB color tuple.
            
        """
        r,g,b,a=self.cm(int(self.m*x))
        # Scale the color values to the range 0-255 and return the color tuple
        return (int(r*self.m),int(g*self.m),int(b*self.m))
 
    def __call__(self, word, font_size, position, orientation,random_state, **kwargs):
        # Returns the color for a given word.
        return self.c(self.f[word])
    

def plot_wordcloud(result_df, score, output_path, figsize=(10, 3), colormap=plt.get_cmap("viridis")):
    label_offset = 0

    df = result_df.copy()

    s = df.loc[:, score].copy().astype(float)
    
    # set the font size based on the ranking of the score
    df["rank"] = s.rank(ascending=asc4rank(score))
    
    # Normalize the score values to be between 0.01 and 0.99
    # if 0 - 1, 0 value will be removed from the word cloud
    df["norm_score"] = preprocessing.minmax_scale(s, feature_range=(0.01,0.99))
    
    # create a color function with the colormap
    cfunc = FreqColorFunc(df.set_index("genename")["norm_score"].to_dict(), colormap)

    # Create a figure with one subplot for the word cloud
    fig = plt.figure(figsize=figsize, dpi=300)
    ax1 = fig.add_subplot(1, 1, 1)
    ax1.axis("off")
    ax1.set_title(f"Rank indes: {score}", fontsize=18)
    
    # Generate the word cloud
    wordcloud = WordCloud(background_color="white", width=800, height=600, max_words=1000).generate_from_frequencies(df.set_index("genename")["rank"].to_dict())
    
    # Recolor the word cloud using the color function
    ax1.imshow(wordcloud.recolor(color_func=cfunc), interpolation="bilinear")
    
    # Add the color bar to the figure
    norm = plt.Normalize(vmin=df.loc[:, score].min(), vmax=df.loc[:, score].max())
    sm = plt.cm.ScalarMappable(cmap=colormap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, orientation='vertical', label=f"{score}", pad=0.01, aspect=15)
    
    # Reduce the number of yticks to only the minimum and maximum
    yticks = [df.loc[:, score].min(), df.loc[:, score].max()]
    yticklabels = ["%.2e" % y for y in yticks]
    cbar.set_ticks(yticks)
    cbar.ax.set_yticklabels(yticklabels)
    
    # Save the figure
    plt.savefig(f"{output_path}/{score}.png", bbox_inches="tight")
    plt.close()
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot wordcloud for DEG')
    parser.add_argument('-i', '--input', help='input file', required=True)
    parser.add_argument('-o', '--output', help='output path', required=True)
    parser.add_argument('-s', '--score', help='score', required=True)
    parser.add_argument("--sep", default="\t", help="A delimiter of input", required=False)
    parser.add_argument('-f', '--figsize', help='figsize', required=False)
    args = parser.parse_args()
    if not os.path.exists(args.output):
        os.makedirs(f"{args.output}/wordcloud", exist_ok=True)
    df = pd.read_csv(args.input, sep=args.sep)
    if args.figsize:
        figsize = tuple([int(i) for i in args.figsize.split(',')])
    else:
        figsize = (10, 3)
    plot_wordcloud(df, args.score, args.output, figsize=figsize)