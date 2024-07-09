import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()
import scipy.sparse as sp
import numpy as np
from .model_MultiGATE import MGATE
from tqdm import tqdm

class MultiGATE():

    def __init__(self, hidden_dims1, hidden_dims2, spot_num, temp, n_epochs=500, lr=0.0001,
                 gradient_clipping=5, nonlinear=True, weight_decay=0.0001,
                 verbose=False, random_seed=2020):
        np.random.seed(random_seed)
        tf.set_random_seed(random_seed)
        self.loss_list = []
        self.lr = lr
        self.n_epochs = n_epochs
        self.gradient_clipping = gradient_clipping
        self.build_placeholders()
        self.verbose = verbose
        
        self.temp = temp
        self.mgate = MGATE(hidden_dims1, hidden_dims2, spot_num, temp, nonlinear, weight_decay)
        self.loss, self.loss_rna, self.loss_atac, self.weight_decay_loss, self.clip_loss, self.H1, self.H2, self.C1, \
        self.C2, self.Cgp, self.ReX1, self.ReX2 = self.mgate(self.A, self.prune_A, self.GP, self.X1, self.X2)
        self.optimize(self.loss)
        self.build_session()
        # self.hidden_dims = hidden_dims

    def build_placeholders(self):
        self.A = tf.sparse_placeholder(dtype=tf.float32)
        self.prune_A = tf.sparse_placeholder(dtype=tf.float32)
        self.GP = tf.sparse_placeholder(dtype=tf.float32)
        self.X1 = tf.placeholder(dtype=tf.float32)
        self.X2 = tf.placeholder(dtype=tf.float32)

    def build_session(self, gpu=True):
        config = tf.ConfigProto()
        config.gpu_options.allow_growth = True
        if gpu == False:
            config.intra_op_parallelism_threads = 0
            config.inter_op_parallelism_threads = 0
        self.session = tf.Session(config=config)
        self.session.run([tf.global_variables_initializer(), tf.local_variables_initializer()])

    def optimize(self, loss):
        optimizer = tf.train.AdamOptimizer(learning_rate=self.lr)
        gradients, variables = zip(*optimizer.compute_gradients(loss))
        gradients, _ = tf.clip_by_global_norm(gradients, self.gradient_clipping)
        self.train_op = optimizer.apply_gradients(zip(gradients, variables))

    def __call__(self, A, prune_A, GP, X1, X2):
        # for epoch in tqdm(range(self.n_epochs)):
        #     self.run_epoch(epoch, A, prune_A, GP, X1, X2)
        
        with tqdm(total=self.n_epochs, desc="Epoch Progress", unit="epoch") as pbar:
            for epoch in range(self.n_epochs):
                loss = self.run_epoch(epoch, A, prune_A, GP, X1, X2)
                pbar.update(1)
                if self.verbose:
                    tqdm.write(f"Epoch: {epoch}, Loss: {loss:.4f}")
                    
    def run_epoch(self, epoch, A, prune_A, GP, X1, X2):

        loss, loss_rna, loss_atac, weight_decay_loss, clip_loss, _ = self.session.run([self.loss, self.loss_rna, self.loss_atac, self.weight_decay_loss, self.clip_loss, self.train_op],
                                         feed_dict={self.A: A,
                                                    self.prune_A: prune_A,
                                                    self.GP: GP,
                                                    self.X1: X1,
                                                    self.X2: X2})
        self.loss_list.append(loss)
        #if self.verbose:
        #    print("Epoch: %s, Loss: %.4f" % (epoch, loss))
        # print("Epoch: %s, Loss: %.4f" % (epoch, loss))
        # print("Epoch: %s, Loss_rna: %.4f" % (epoch, loss_rna))
        # print("Epoch: %s, Loss_atac: %.4f" % (epoch, loss_atac))
        # print("Epoch: %s, Loss_weight_decay: %.4f" % (epoch, weight_decay_loss))
        # print("Epoch: %s, Loss_clip: %.4f" % (epoch, clip_loss))
        return loss

    def infer(self, A, prune_A, GP, X1, X2):
        H1, H2, C1, C2, Cgp, ReX1, ReX2 = self.session.run([self.H1, self.H2, self.C1, self.C2, self.Cgp, self.ReX1, self.ReX2],
                           feed_dict={self.A: A,
                                      self.prune_A: prune_A,
                                      self.GP: GP,
                                      self.X1: X1,
                                      self.X2: X2})

        return H1, H2, self.Conbine_Atten_l(C1), self.Conbine_Atten_l(C2), self.Conbine_Atten_l(Cgp), self.loss_list, ReX1, ReX2

    def Conbine_Atten_l(self, input):
        return [sp.coo_matrix((input[layer][1], (input[layer][0][:, 0], input[layer][0][:, 1])), shape=(input[layer][2][0], input[layer][2][1])) for layer in input]
   