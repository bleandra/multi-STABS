###############################################################################
# RESULTS PLOTTING
#
# this script counts the individuals in each state per timestep
# and visualises the results of the model run
###############################################################################



# arranging data---------------------------------------------------------------


# total malaria
ordered_human_stack<-human_stack[order(human_stack$t),]
# because of the multispecies infections, one individual can be in two states at once (one for each species)

Eh<-cumsum((1*(ordered_human_stack$Pf_state_current==1) + 1*(ordered_human_stack$Pv_state_current==1) -
            (1*(ordered_human_stack$Pf_state_prev==1) + 1*(ordered_human_stack$Pv_state_prev==1) )))
Ih<-initial_I+cumsum(1*(ordered_human_stack$Pf_state_current==2) + 1*(ordered_human_stack$Pv_state_current==2) - 
                     (1*(ordered_human_stack$Pf_state_prev==2) + 1*(ordered_human_stack$Pv_state_prev==2) ))
Gh<-cumsum(1*(ordered_human_stack$Pf_state_current==4) + 1*(ordered_human_stack$Pv_state_current==4) -
           (1*(ordered_human_stack$Pf_state_prev==4) + 1*(ordered_human_stack$Pv_state_prev==4)))
Th<-cumsum(1*(ordered_human_stack$Pf_state_current==5) + 1*(ordered_human_stack$Pv_state_current==5) -
           (1*(ordered_human_stack$Pf_state_prev==5) + 1*(ordered_human_stack$Pv_state_prev==5) ))
Lh<-cumsum(1*(ordered_human_stack$Pf_state_current==6) + 1*(ordered_human_stack$Pv_state_current==6) -
          ( 1*(ordered_human_stack$Pf_state_prev==6) + 1*(ordered_human_stack$Pv_state_prev==6) ))
Sh<- (tot_h - Eh - Ih - Gh -Th -Lh)
human_I<-data.frame(cbind(ordered_human_stack$t[1:h_stack_count],
                          Sh[1:h_stack_count],
                          Eh[1:h_stack_count],
                          Ih[1:h_stack_count], 
                          Gh[1:h_stack_count],
                          Th[1:h_stack_count],
                          Lh[1:h_stack_count]))
names(human_I)<-c("t","Sh","Eh","Ih","Gh", "Th", "Lh")
human_I<-human_I[which(human_I$t<max_t),]



# plots -----------------------------------------------------
ggplot(human_I)+geom_step(aes(t,Sh))


ordered_m_stack<-mosquito_stack[order(mosquito_stack$t),]
ordered_m_stack<-ordered_m_stack[1:stack_size,]

state_vec<-data.frame(matrix(NA,ncol=6,nrow=3*m_count))
names(state_vec)<-c('tm','event','Em','Im','Rm', 'species')
row_state<-0;
for(ii in 1:m_count){
  #data relating to each mosquito
  mozzy_ii<-which(ordered_m_stack$mosquito_id==ii)
  
    #recording the exposure
  state_vec[row_state+1,1]<-ordered_m_stack$t[mozzy_ii[1]]
  state_vec[row_state+1,2]<-1
  state_vec[row_state+1,6]<-ordered_m_stack$mal_species[mozzy_ii[1]]
  
  #recording the removal
  tG_index<-which(ordered_m_stack$event[mozzy_ii]==3)
  if(length(tG_index)==0){
    print(ii)
    state_vec[row_state+3,1]<-max_t
  }else{
    state_vec[row_state+3,1]<-ordered_m_stack$t[mozzy_ii[tG_index]]
  }
  state_vec[row_state+3,2]<-3
  state_vec[row_state+3,6]<-ordered_m_stack$mal_species[mozzy_ii[1]]
  
  #recording infectiousness
  state_vec[row_state+2,1]<-min(state_vec[row_state+3,1],exposure_time[ii])
  #note the infectiousness time effectively happens at the same time as the recovery time if removal occures prior to infectiousness
  state_vec[row_state+2,2]<-2
  state_vec[row_state+2,6]<-ordered_m_stack$mal_species[mozzy_ii[1]]
  
  row_state<-row_state+3
}
state_vec<-state_vec[order(state_vec[,1]),]
state_vec_pf<-state_vec
state_vec_pv<-state_vec
state_vec$Em<-cumsum(1*(state_vec$event==1)-1*(state_vec$event==2))
state_vec$Im<-cumsum(1*(state_vec$event==2)-1*(state_vec$event==3))
#R is irrelevant in infinite population of mosquitos
state_vec$Rm<-cumsum(1*(state_vec$event==3))

#state_vec<-state_vec[1:row_state,]
state_vec<-state_vec[which(state_vec$tm<=max_t-1),] # bc we're not using the removal of the mozzies in the plots this is sufficient


stuff<-data.frame(matrix(NA,ncol=length(human_I[1,]),nrow=(length(state_vec[,1])-length(human_I[,1]))))
names(stuff)<-names(human_I)
big_df<-cbind(state_vec,rbind(human_I,stuff))

g1<-ggplot(big_df)+
  geom_step(aes(tm,Em,colour="Mosquito",linetype="Exposed",alpha="Exposed",size="Exposed"))+
  geom_step(aes(tm,Im,colour="Mosquito",linetype="Infectious",alpha="Infectious",size="Infectious"))+
  geom_step(aes(t,Eh,colour="Human",linetype="Exposed",alpha="Exposed",size="Exposed"))+
  geom_step(aes(t,Ih,colour="Human",linetype="Infectious",alpha="Infectious",size="Infectious"))+
  geom_step(aes(t,Th,colour="Human",linetype="Standard Treatment",alpha="Standard Treatment",size="Standard Treatment"))+
  geom_step(aes(t,Gh,colour="Human",linetype="Radical Cure",alpha="Radical Cure",size="Radical Cure"))+
  geom_step(aes(t,Lh,colour="Human",linetype="Latent Hypnozoites",alpha="Latent Hypnozoites",size="Latent Hypnozoites"))+
  labs(x = "time (days)",y= "Number of Individuals")+
  scale_linetype_manual(values = c('Exposed' = "dashed",'Infectious' = "solid", 'Standard Treatment' = "dotted", 'Radical Cure' = "twodash", 'Latent Hypnozoites' = "dotdash")) +
  scale_alpha_manual(values = c('Exposed' = 0.5,'Infectious' = 1, 'Standard Treatment' = 0.75, 'Radical Cure' = 0.75, 'Latent Hypnozoites' = 0.75)) +
  scale_size_manual(values = c('Exposed' = 0.5,'Infectious' = 1.1, 'Standard Treatment' = 0.5, 'Radical Cure' = 0.5, 'Latent Hypnozoites' = 0.5)) +
  guides(alpha = "none",size="none") +
  ggtitle("all Malaria")

g1
# ggsave("output/multi_STABS_HM.png",g1, width = 5, height = 5)


comp_plot<-ggplot(big_df)+
  geom_step(aes(t,Eh,colour="Exposed",alpha="Exposed",size="Exposed"))+
  geom_step(aes(t,Ih,colour="Infectious",alpha="Infectious",size="Infectious"))+
  geom_step(aes(t,Th,colour="Standard Treatment",alpha="Standard Treatment",size="Standard Treatment"))+
  geom_step(aes(t,Gh,colour="Radical Cure",alpha="Radical Cure",size="Radical Cure"))+
  geom_step(aes(t,Lh,colour="Latent Hypnozoites",alpha="Latent Hypnozoites",size="Latent Hypnozoites"))+
  labs(x = "time (days)",y= "Number of Individuals")+
  #scale_colour_manual(values = c('Mosquito' = "grey66",'Human' = "red3")) +
  # scale_linetype_manual(values = c('Exposed' = "dashed",'Infectious' = "solid", 'Standard Treatment' = "dotted", 'Radical Cure' = "twodash", 'Latent Hypnozoites' = "dotdash")) +
  scale_alpha_manual(values = c('Exposed' = 0.5,'Infectious' = 1, 'Standard Treatment' = 0.75, 'Radical Cure' = 0.75, 'Latent Hypnozoites' = 0.75)) +
  scale_size_manual(values = c('Exposed' = 0.5,'Infectious' = 1.1, 'Standard Treatment' = 0.5, 'Radical Cure' = 0.5, 'Latent Hypnozoites' = 0.5)) +
  guides(alpha = "none",size="none") +
  ggtitle("Human States For All Malaria")

comp_plot
ggsave("output/multiSTABS_human_comps.png",plot = comp_plot, width = 5, height = 5)


# species-specific -----------------------------------------------------------------------------
# 
# # for pf#--------------------------------------------------------------------
Eh_pf<-cumsum((1*(ordered_human_stack$Pf_state_current==1)  -
              (1*(ordered_human_stack$Pf_state_prev==1)  )))
Ih_pf<-initial_I_pf+cumsum(1*(ordered_human_stack$Pf_state_current==2) - 
                       (1*(ordered_human_stack$Pf_state_prev==2) ))
Gh_pf<-cumsum(1*(ordered_human_stack$Pf_state_current==4) -
             (1*(ordered_human_stack$Pf_state_prev==4) ))
Th_pf<-cumsum(1*(ordered_human_stack$Pf_state_current==5) -
             (1*(ordered_human_stack$Pf_state_prev==5) ))
Lh_pf<-cumsum(1*(ordered_human_stack$Pf_state_current==6) - # can be taken out as it should be 0 
             ( 1*(ordered_human_stack$Pf_state_prev==6) ))
Sh_pf<- (tot_h - Eh_pf - Ih_pf - Gh_pf -Th_pf -Lh_pf) # this is an approximation
human_I_pf<-data.frame(cbind(ordered_human_stack$t[1:h_stack_count],
                          Sh_pf[1:h_stack_count],
                          Eh_pf[1:h_stack_count],
                          Ih_pf[1:h_stack_count], 
                          Gh_pf[1:h_stack_count],
                          Th_pf[1:h_stack_count],
                          Lh_pf[1:h_stack_count])) # this compartment can be taken out
names(human_I_pf)<-c("t","Sh_pf","Eh_pf","Ih_pf","Gh_pf", "Th_pf", "Lh_pf")
human_I_pf<-human_I_pf[which(human_I_pf$t<max_t),]

names(state_vec_pf)<-c('tm','event','Em_pf','Im_pf','Rm_pf', 'species')
state_vec_pf$Em_pf<-cumsum(1*(state_vec_pf$event==1& state_vec_pf$species==1)-
                             1*(state_vec_pf$event==2& state_vec_pf$species==1))
state_vec_pf$Im_pf<-cumsum(1*(state_vec_pf$event==2& state_vec_pf$species==1)-
                             1*(state_vec_pf$event==3& state_vec_pf$species==1))
#R is irrelevant in infinite population of mosquitos
state_vec_pf$Rm_pf<-cumsum(1*(state_vec_pf$event==3& state_vec_pf$species==1))

state_vec_pf<-state_vec_pf[which(state_vec_pf$tm<=max_t-1),] # bc we're not using the removal of the mozzies in the plots this is sufficient


stuff_pf<-data.frame(matrix(NA,ncol=length(human_I_pf[1,]),nrow=(length(state_vec_pf[,1])-length(human_I_pf[,1]))))
names(stuff_pf)<-names(human_I_pf)
big_df_pf<-cbind(state_vec_pf,rbind(human_I_pf,stuff_pf))


g1_pf<-ggplot(big_df_pf)+
  geom_step(aes(tm,Em_pf,colour="Mosquito",linetype="Exposed",alpha="Exposed",size="Exposed"))+
  geom_step(aes(tm,Im_pf,colour="Mosquito",linetype="Infectious",alpha="Infectious",size="Infectious"))+
  geom_step(aes(t,Eh_pf,colour="Human",linetype="Exposed",alpha="Exposed",size="Exposed"))+
  geom_step(aes(t,Ih_pf,colour="Human",linetype="Infectious",alpha="Infectious",size="Infectious"))+
  geom_step(aes(t,Th_pf,colour="Human",linetype="Standard Treatment",alpha="Standard Treatment",size="Standard Treatment"))+
  geom_step(aes(t,Gh_pf,colour="Human",linetype="Radical Cure",alpha="Radical Cure",size="Radical Cure"))+
  # geom_step(aes(t,Lh_pf,colour="Human",linetype="Latent Hypnozoites",alpha="Latent Hypnozoites",size="Latent Hypnozoites"))+
  labs(x = "time (days)",y= "Number of Individuals")+
  #scale_colour_manual(values = c('Mosquito' = "grey66",'Human' = "red3")) +
  scale_linetype_manual(values = c('Exposed' = "dashed",'Infectious' = "solid", 'Standard Treatment' = "dotted", 'Radical Cure' = "twodash")) +
  scale_alpha_manual(values = c('Exposed' = 0.5,'Infectious' = 1, 'Standard Treatment' = 0.75, 'Radical Cure' = 0.75)) +
  scale_size_manual(values = c('Exposed' = 0.5,'Infectious' = 1.1, 'Standard Treatment' = 0.5, 'Radical Cure' = 0.5)) +
  guides(alpha = "none",size="none") +
  ggtitle("P. falciparum")

g1_pf
# ggsave("output/multi_STABS_HM_pf.png",g1_pf, width = 5, height = 5)




# # for pv #--------------------------------------------------------------------
Eh_pv<-cumsum((1*(ordered_human_stack$Pv_state_current==1)  -
                 (1*(ordered_human_stack$Pv_state_prev==1)  )))
Ih_pv<-initial_I_pv+cumsum(1*(ordered_human_stack$Pv_state_current==2) - 
                             (1*(ordered_human_stack$Pv_state_prev==2) ))
Gh_pv<-cumsum(1*(ordered_human_stack$Pv_state_current==4) -
                (1*(ordered_human_stack$Pv_state_prev==4) ))
Th_pv<-cumsum(1*(ordered_human_stack$Pv_state_current==5) -
                (1*(ordered_human_stack$Pv_state_prev==5) ))
Lh_pv<-cumsum(1*(ordered_human_stack$Pv_state_current==6) -
                ( 1*(ordered_human_stack$Pv_state_prev==6) ))
Sh_pv<- (tot_h - Eh_pv - Ih_pv - Gh_pv -Th_pv -Lh_pv) # this is an approximation
human_I_pv<-data.frame(cbind(ordered_human_stack$t[1:h_stack_count],
                             Sh_pv[1:h_stack_count],
                             Eh_pv[1:h_stack_count],
                             Ih_pv[1:h_stack_count], 
                             Gh_pv[1:h_stack_count],
                             Th_pv[1:h_stack_count],
                             Lh_pv[1:h_stack_count]))
names(human_I_pv)<-c("t","Sh_pv","Eh_pv","Ih_pv","Gh_pv", "Th_pv", "Lh_pv")
human_I_pv<-human_I_pv[which(human_I_pv$t<max_t),]

names(state_vec_pv)<-c('tm','event','Em_pv','Im_pv','Rm_pv', 'species')
state_vec_pv$Em_pv<-cumsum(1*(state_vec_pv$event==1& state_vec_pv$species==3)-
                             1*(state_vec_pv$event==2& state_vec_pv$species==3))
state_vec_pv$Im_pv<-cumsum(1*(state_vec_pv$event==2& state_vec_pv$species==3)-
                             1*(state_vec_pv$event==3& state_vec_pv$species==3))
#R is irrelevant in infinite population of mosquitos
state_vec_pv$Rm_pv<-cumsum(1*(state_vec_pv$event==3& state_vec_pv$species==3))

state_vec_pv<-state_vec_pv[which(state_vec_pv$tm<=max_t-1),] # bc we're not using the removal of the mozzies in the plots this is sufficient


stuff_pv<-data.frame(matrix(NA,ncol=length(human_I_pv[1,]),nrow=(length(state_vec_pv[,1])-length(human_I_pv[,1]))))
names(stuff_pv)<-names(human_I_pv)
big_df_pv<-cbind(state_vec_pv,rbind(human_I_pv,stuff_pv))

g1_pv<-ggplot(big_df_pv)+
  geom_step(aes(tm,Em_pv,colour="Mosquito",linetype="Exposed",alpha="Exposed",size="Exposed"))+
  geom_step(aes(tm,Im_pv,colour="Mosquito",linetype="Infectious",alpha="Infectious",size="Infectious"))+
  geom_step(aes(t,Eh_pv,colour="Human",linetype="Exposed",alpha="Exposed",size="Exposed"))+
  geom_step(aes(t,Ih_pv,colour="Human",linetype="Infectious",alpha="Infectious",size="Infectious"))+
  geom_step(aes(t,Th_pv,colour="Human",linetype="Standard Treatment",alpha="Standard Treatment",size="Standard Treatment"))+
  geom_step(aes(t,Gh_pv,colour="Human",linetype="Radical Cure",alpha="Radical Cure",size="Radical Cure"))+
  geom_step(aes(t,Lh_pv,colour="Human",linetype="Latent Hypnozoites",alpha="Latent Hypnozoites",size="Latent Hypnozoites"))+
  labs(x = "time (days)",y= "Number of Individuals")+
  #scale_colour_manual(values = c('Mosquito' = "grey66",'Human' = "red3")) +
  scale_linetype_manual(values = c('Exposed' = "dashed",'Infectious' = "solid", 'Standard Treatment' = "dotted", 'Radical Cure' = "twodash", 'Latent Hypnozoites' = "dotdash")) +
  scale_alpha_manual(values = c('Exposed' = 0.5,'Infectious' = 1, 'Standard Treatment' = 0.75, 'Radical Cure' = 0.75, 'Latent Hypnozoites' = 0.75)) +
  scale_size_manual(values = c('Exposed' = 0.5,'Infectious' = 1.1, 'Standard Treatment' = 0.5, 'Radical Cure' = 0.5, 'Latent Hypnozoites' = 0.5)) +
  guides(alpha = "none",size="none") +
  ggtitle("P. vivax")

g1_pv
# ggsave("output/multi_STABS_HM_pv.png",plot=g1_pv, width = 5, height = 5)


## get plots arranged #---------------------------------------------------------

(g1 + theme(legend.position="none")) + (g1_pf+ theme(legend.position="none")) + g1_pv + 
  plot_annotation(title = current_scen) & 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("output/multi_STABS_grid.png",width = 15, height = 5)
